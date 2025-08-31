'use client';

import React, { useCallback, useEffect, useMemo, useState } from 'react';
import axios from 'axios';
import { Card, CardBody, CardFooter, CardHeader } from '@nextui-org/card';
import { Button } from '@nextui-org/button';
import { Progress } from '@nextui-org/progress';
import { Chip } from '@nextui-org/chip';
import { Tooltip } from '@nextui-org/tooltip';
import { Divider } from '@nextui-org/divider';
import { Slider } from '@nextui-org/slider';
import { Modal, ModalBody, ModalContent, ModalFooter } from '@nextui-org/modal';
import { Spinner } from '@nextui-org/spinner';

const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

type SliceItem = { slice_index: number; image: string };
type InlineResult = {
  segmented_slices: SliceItem[];
  download_url: string;
  volumes?: { pancreas_voxels: number; lesion_voxels: number; total_seg_voxels: number };
  status?: string;
};
type JobSubmit = { job_id: string; status_url: string; download_url: string };
type JobStatus = { status: 'queued' | 'running' | 'done' | 'error'; error?: string; mask_path?: string };

const WINDOW_PRESETS: { label: string; wl: number; ww: number }[] = [
  // presets commonly used in abdomen CT (approximate visual mapping)
  { label: 'Soft Tissue', wl: 50, ww: 350 },
  { label: 'Liver', wl: 60, ww: 150 },
  { label: 'Pancreas', wl: 60, ww: 300 },
  { label: 'Bone', wl: 300, ww: 1500 },
];

export default function PancreasCTSegmentation() {
  // inputs
  const [file, setFile] = useState<File | null>(null);
  const isFormValid = !!file;

  // run state
  const [loading, setLoading] = useState(false);
  const [progress, setProgress] = useState(0);
  const [error, setError] = useState<string | null>(null);
  const [lastRunAt, setLastRunAt] = useState<string | null>(null);

  // results (inline)
  const [result, setResult] = useState<InlineResult | null>(null);
  const slices = result?.segmented_slices ?? [];
  const sliceCount = slices.length;
  const [idx, setIdx] = useState(0);

  // optional job mode
  const [jobId, setJobId] = useState<string | null>(null);
  const [jobDownloadUrl, setJobDownloadUrl] = useState<string | null>(null);

  // viewer controls (simple brightness/contrast; WL/WW presets map to these)
  const [brightness, setBrightness] = useState(100);
  const [contrast, setContrast] = useState(100);
  const [zoomOpen, setZoomOpen] = useState(false);

  // health (optional; skip if you haven’t added a health endpoint)
  const [healthChecking, setHealthChecking] = useState(false);
  const [healthOk, setHealthOk] = useState<boolean | null>(null);

  useEffect(() => {
    // no /health route in your pancreas router yet — leave off by default
    if (!process.env.NEXT_PUBLIC_PANCREAS_HEALTH_URL) return;
    const url = process.env.NEXT_PUBLIC_PANCREAS_HEALTH_URL!;
    setHealthChecking(true);
    axios
      .get(url, { timeout: 5000 })
      .then(() => setHealthOk(true))
      .catch(() => setHealthOk(false))
      .finally(() => setHealthChecking(false));
  }, []);

  // smooth fake progress (UX)
  useEffect(() => {
    if (!loading) return;
    let raf = 0;
    const start = performance.now();
    const duration = 120_000; // 2 min visual ramp; backend may take longer
    const tick = (t: number) => {
      const x = Math.min((t - start) / duration, 1);
      const eased = 1 - (1 - x) * (1 - x);
      setProgress(Math.min(95, Math.round(eased * 100)));
      if (x < 1) raf = requestAnimationFrame(tick);
    };
    raf = requestAnimationFrame(tick);
    return () => cancelAnimationFrame(raf);
  }, [loading]);

  const formatBytes = (bytes?: number) => {
    if (!bytes && bytes !== 0) return '';
    const k = 1024, sizes = ['B', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return `${(bytes / Math.pow(k, i)).toFixed(2)} ${sizes[i]}`;
  };
  const pct = (p: number, t: number) => (t ? (p / t) * 100 : 0);

  const applyPreset = useCallback((wl: number, ww: number) => {
    // crude mapping WL/WW → CSS brightness/contrast to mimic windowing on a baked PNG
    // (not physically accurate; just a visual aid)
    const c = Math.min(200, Math.max(50, Math.round((ww / 350) * 100))); // reference WW≈350
    const b = Math.min(150, Math.max(50, Math.round(100 + (wl - 50) * 0.5)));
    setBrightness(b);
    setContrast(c);
  }, []);

  const startSegmentation = useCallback(async () => {
    if (!isFormValid || !file) {
      setError('Please upload a CT volume (NIfTI .nii or .nii.gz).');
      return;
    }
    setError(null);
    setLoading(true);
    setProgress(0);
    setResult(null);
    setIdx(0);
    setJobId(null);
    setJobDownloadUrl(null);

    const formData = new FormData();
    formData.append('file', file);

    try {
      const r = await axios.post(`${API_BASE_URL}/imaging/pancreas-tumor/segment`, formData, {
        headers: { 'Content-Type': 'multipart/form-data' },
        timeout: 600_000, // 10 minutes
      });

      // Handle either inline or job-style response
      if (r.data?.segmented_slices) {
        const data = r.data as InlineResult;
        setResult(data);
        setProgress(100);
        setLastRunAt(new Date().toISOString());
        if (data.segmented_slices.length) setIdx(0);
      } else if (r.data?.job_id) {
        const data = r.data as JobSubmit;
        setJobId(data.job_id);
        setJobDownloadUrl(data.download_url);
        // poll status
        await pollJob(data.status_url);
      } else {
        throw new Error('Unexpected API response format.');
      }
    } catch (err: any) {
      const detail = err?.response?.data?.detail || err?.message || 'Unknown error';
      setError(`Segmentation failed: ${detail}`);
    } finally {
      setLoading(false);
    }
  }, [file, isFormValid]);

  const pollJob = async (statusUrl: string) => {
    let tries = 0;
    const maxTries = 600; // ~30–50 min if every 3–5s (tune as needed)
    while (tries < maxTries) {
      tries += 1;
      try {
        const s = await axios.get<JobStatus>(`${API_BASE_URL}${statusUrl}`, { timeout: 15_000 });
        if (s.data.status === 'done') {
          setProgress(100);
          setLastRunAt(new Date().toISOString());
          // job mode doesn’t include overlays; offer download when ready
          return;
        }
        if (s.data.status === 'error') {
          setError(s.data.error || 'Job failed.');
          return;
        }
      } catch (e) {
        // ignore transient errors
      }
      await new Promise((res) => setTimeout(res, 4000));
    }
    setError('Timed out waiting for job to complete.');
  };

  const pancreasVox = result?.volumes?.pancreas_voxels ?? 0;
  const lesionVox = result?.volumes?.lesion_voxels ?? 0;
  const totalVox = result?.volumes?.total_seg_voxels ?? (pancreasVox + lesionVox);

  return (
    <div className="max-w-7xl mx-auto p-4 md:p-8 bg-gray-50 min-h-screen">
      {/* Header */}
      <div className="sticky top-0 z-10 mb-6">
        <Card className="border border-gray-200 shadow-sm">
          <CardBody className="py-5">
            <div className="flex flex-col gap-6 md:flex-row md:items-center md:justify-between">
              <div>
                <h1 className="text-2xl md:text-3xl font-semibold tracking-tight text-gray-900">
                  Pancreas CT Segmentation
                </h1>
                <p className="text-sm text-gray-600 mt-1">Upload a single CT volume (NIfTI) to segment pancreas ± lesion.</p>
                <div className="mt-2 flex flex-wrap gap-2">
                  <Chip variant="flat" color="primary">DiNTS</Chip>
                  <Chip variant="flat">NIfTI .nii / .nii.gz</Chip>
                </div>
              </div>
              <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                <Kpi label="Service" value={healthChecking ? 'Checking…' : healthOk === false ? 'Offline' : 'Online'} tone={healthOk === false ? 'bad' : 'ok'} />
                <Kpi label="File" value={file ? 'Loaded' : '—'} />
                <Kpi label="Slices" value={sliceCount ? String(sliceCount) : jobId ? '—' : '—'} />
                <Kpi label="Last run" value={lastRunAt ? new Date(lastRunAt).toLocaleString() : '—'} />
              </div>
            </div>
          </CardBody>
        </Card>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-12 gap-6">
        {/* Left: input & progress */}
        <section className="lg:col-span-3">
          <Card className="border border-gray-200 shadow-sm">
            <CardHeader className="bg-gray-50"><p className="text-lg font-semibold">1 · Input CT Volume</p></CardHeader>
            <CardBody>
              <Tooltip content="NIfTI (.nii / .nii.gz) CT volume" placement="top">
                <div className={`rounded-md border ${file ? 'border-emerald-300' : 'border-gray-200'} bg-white p-3`}>
                  <div className="flex items-center justify-between mb-2">
                    <span className="font-medium">CT NIfTI</span>
                    {file ? <Chip size="sm" color="success" variant="flat">Loaded</Chip> : <Chip size="sm" variant="flat">Awaiting</Chip>}
                  </div>
                  <input
                    type="file"
                    accept=".nii,.nii.gz"
                    className="block w-full text-sm text-gray-600 file:mr-2 file:py-1 file:px-3 file:rounded-md file:border-0 file:bg-indigo-50 file:text-indigo-700 hover:file:bg-indigo-100"
                    onChange={(e) => setFile(e.target.files?.[0] ?? null)}
                    disabled={loading}
                  />
                  {file && (
                    <div className="mt-1 text-[11px] text-gray-500 break-all">
                      {file.name} • {formatBytes(file.size)}
                    </div>
                  )}
                </div>
              </Tooltip>

              <div className="mt-5 flex flex-col items-stretch gap-2">
                <Button
                  color="primary"
                  isDisabled={!isFormValid || loading}
                  isLoading={loading}
                  onPress={startSegmentation}
                >
                  {loading ? 'Processing…' : 'Run Segmentation'}
                </Button>
                {error && <div className="rounded-md border-l-4 border-red-600 bg-red-50 p-3 text-sm text-red-800">{error}</div>}
                {jobId && !result && (
                  <div className="rounded-md border border-amber-200 bg-amber-50 p-3 text-xs text-amber-800">
                    Job submitted: <span className="font-mono">{jobId}</span>. This may take a few minutes.
                  </div>
                )}
              </div>
            </CardBody>
          </Card>

          <Card className="mt-6 border border-gray-200 shadow-sm">
            <CardHeader className="bg-gray-50"><p className="text-lg font-semibold">2 · Processing Status</p></CardHeader>
            <CardBody>
              {loading ? (
                <>
                  <div className="flex justify-between text-sm mb-1">
                    <span>Segmentation in progress</span><span>{progress}%</span>
                  </div>
                  <Progress value={progress} color="primary" className="h-2" />
                  <div className="text-xs text-gray-600 mt-3">Please keep this tab open.</div>
                </>
              ) : (
                <p className="text-sm text-gray-500">{result ? 'Completed.' : 'No analysis in progress.'}</p>
              )}
            </CardBody>
          </Card>
        </section>

        {/* Center: viewer */}
        <section className="lg:col-span-6">
          <Card className="border border-gray-200 shadow-sm">
            <CardHeader className="bg-gray-50 flex items-center justify-between">
              <p className="font-semibold">3 · Viewer</p>
              <div className="flex items-center gap-3 text-xs text-gray-600">
                <span>WL/WW</span>
                <div className="hidden md:flex gap-2">
                  {WINDOW_PRESETS.map((p) => (
                    <Button key={p.label} size="sm" variant="flat" onPress={() => applyPreset(p.wl, p.ww)}>
                      {p.label}
                    </Button>
                  ))}
                </div>
                <Button size="sm" variant="flat" onPress={() => setZoomOpen(true)}>Open Lightbox</Button>
              </div>
            </CardHeader>
            <CardBody>
              {!result && !jobId ? (
                <div className="h-[440px] flex items-center justify-center text-gray-500">Upload and run to view slices.</div>
              ) : jobId && !result ? (
                <div className="h-[440px] flex flex-col items-center justify-center gap-3 text-gray-600">
                  <Spinner color="primary" />
                  <div className="text-center text-sm">
                    Job in progress. You’ll get a download link when ready.
                  </div>
                  {jobDownloadUrl && (
                    <Button as="a" href={jobDownloadUrl} variant="bordered" color="primary">
                      Download (when ready)
                    </Button>
                  )}
                </div>
              ) : (
                <>
                  <div className="relative w-full h-[440px] rounded-md overflow-hidden bg-black">
                    <img
                      src={slices[idx]?.image}
                      alt={`Slice ${slices[idx]?.slice_index}`}
                      className="h-full w-full object-contain"
                      style={{ filter: `brightness(${brightness}%) contrast(${contrast}%)` }}
                    />
                  </div>

                  <div className="mt-4 grid grid-cols-1 md:grid-cols-2 gap-4">
                    <div>
                      <div className="flex items-center justify-between text-xs text-gray-600 mb-1">
                        <span>Slice</span>
                        <span>#{slices[idx]?.slice_index ?? '—'}</span>
                      </div>
                      <Slider
                        aria-label="Slice selector"
                        minValue={0}
                        maxValue={Math.max(0, sliceCount - 1)}
                        value={idx}
                        onChange={(v) => setIdx(Number(v))}
                        size="sm"
                      />
                    </div>
                    <div>
                      <div className="flex items-center justify-between text-xs text-gray-600 mb-1">
                        <span>Brightness / Contrast</span>
                        <span>{brightness}% / {contrast}%</span>
                      </div>
                      <div className="grid grid-cols-2 gap-3">
                        <Slider minValue={50} maxValue={150} value={brightness} onChange={(v)=>setBrightness(Number(v))} aria-label="Brightness" size="sm" />
                        <Slider minValue={50} maxValue={200} value={contrast}  onChange={(v)=>setContrast(Number(v))}  aria-label="Contrast"  size="sm" />
                      </div>
                    </div>
                  </div>

                  <Divider className="my-4" />

                  <div className="flex flex-wrap items-center gap-4 text-xs">
                    <Legend swatch="bg-green-500" label="Pancreas" />
                    <Legend swatch="bg-red-500"   label="Lesion" />
                  </div>

                  {/* Quick jump */}
                  {sliceCount > 0 && (
                    <div className="mt-4">
                      <p className="text-sm font-medium mb-2">Quick Slice Jump</p>
                      <div className="flex gap-2 overflow-x-auto pb-1">
                        {slices.map((s, i) => (
                          <Button
                            key={s.slice_index}
                            size="sm"
                            variant={i === idx ? 'solid' : 'flat'}
                            color={i === idx ? 'primary' : 'default'}
                            onPress={() => setIdx(i)}
                          >
                            #{s.slice_index}
                          </Button>
                        ))}
                      </div>
                    </div>
                  )}
                </>
              )}
            </CardBody>
          </Card>
        </section>

        {/* Right: metrics & download */}
        <section className="lg:col-span-3">
          <Card className="border border-gray-200 shadow-sm">
            <CardHeader className="bg-gray-50"><p className="font-semibold">4 · Metrics</p></CardHeader>
            <CardBody>
              {!result ? (
                <p className="text-sm text-gray-500">
                  Metrics will appear after segmentation.
                  {jobId && !result ? ' (Job mode: only NIfTI download available.)' : ''}
                </p>
              ) : (
                <>
                  <div className="grid grid-cols-2 gap-3 mb-4">
                    <Metric label="Total voxels" value={fmtVox(totalVox)} />
                    <Metric label="Pancreas voxels" value={fmtVox(pancreasVox)} tone="green" />
                    <Metric label="Lesion voxels" value={fmtVox(lesionVox)} tone="red" />
                    <Metric label="Lesion %" value={`${Math.round(pct(lesionVox, totalVox))}%`} tone="red" />
                  </div>
                </>
              )}
            </CardBody>
            {(result?.download_url || jobDownloadUrl) && (
              <CardFooter className="bg-gray-50">
                <Button
                  as="a"
                  href={result?.download_url ?? jobDownloadUrl ?? '#'}
                  color="primary"
                  variant="bordered"
                  className="w-full"
                >
                  Download Segmentation (NIfTI)
                </Button>
              </CardFooter>
            )}
          </Card>
        </section>
      </div>

      {/* Lightbox */}
      <Modal isOpen={zoomOpen} onClose={() => setZoomOpen(false)} size="5xl">
        <ModalContent>
          <ModalBody className="p-0 bg-black">
            {slices[idx] && (
              <img
                src={slices[idx].image}
                alt="Zoomed slice"
                className="w-full h-full object-contain"
                style={{ filter: `brightness(${brightness}%) contrast(${contrast}%)` }}
              />
            )}
          </ModalBody>
          <ModalFooter className="justify-end">
            <Button onPress={() => setZoomOpen(false)}>Close</Button>
          </ModalFooter>
        </ModalContent>
      </Modal>
    </div>
  );
}

/* ---------- helpers ---------- */
function Kpi({ label, value, tone }: { label: string; value: string; tone?: 'ok' | 'bad' | 'base' }) {
  const color = tone === 'ok' ? 'text-emerald-700' : tone === 'bad' ? 'text-red-700' : 'text-gray-900';
  return (
    <div className="rounded-lg border border-gray-200 bg-white p-3 text-center">
      <p className="text-[11px] uppercase tracking-wide text-gray-500">{label}</p>
      <p className={`font-mono text-xl font-semibold ${color}`}>{value}</p>
    </div>
  );
}
function Metric({ label, value, tone = 'base' }: { label: string; value: string; tone?: 'base' | 'red' | 'green' }) {
  const map: Record<string, string> = {
    base: 'text-gray-900',
    red: 'text-red-600',
    green: 'text-green-600',
  };
  return (
    <div className="rounded-lg bg-gray-50 p-3 text-center">
      <div className={`text-xl font-semibold ${map[tone]}`}>{value}</div>
      <div className="text-xs text-gray-600">{label}</div>
    </div>
  );
}
function Legend({ swatch, label }: { swatch: string; label: string }) {
  return (
    <div className="flex items-center gap-2">
      <span className={`h-3 w-3 rounded-full ${swatch}`} />
      <span className="text-xs text-gray-700">{label}</span>
    </div>
  );
}
function fmtVox(v: number) {
  return v < 1000 ? `${v}` : `${(v / 1000).toFixed(1)}k`;
}
