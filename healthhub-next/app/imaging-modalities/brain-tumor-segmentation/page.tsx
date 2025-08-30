'use client';

import React, { useState, useCallback, useEffect, useMemo } from 'react';
import axios from 'axios';
import { Card, CardBody, CardFooter, CardHeader } from '@nextui-org/card';
import { Button } from '@nextui-org/button';
import { Progress } from '@nextui-org/progress';
import { Tooltip } from '@nextui-org/tooltip';
import { Chip } from '@nextui-org/chip';
import { Modal, ModalContent, ModalBody, ModalFooter } from '@nextui-org/modal';
import { Spinner } from '@nextui-org/spinner';
import { Slider } from '@nextui-org/slider';
import { Divider } from '@nextui-org/divider';
import { Tabs, Tab } from '@nextui-org/tabs';

const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

type ModalityType = 'T1' | 'T1c' | 'T2' | 'FLAIR';
type Modalities = Record<ModalityType, File | null>;

interface SegmentationSlice {
  slice_index: number;
  image: string; // server returns baked overlay image
}

interface SegmentationResult {
  segmented_slices: SegmentationSlice[];
  download_url: string;
  tumor_volumes: { necrotic: number; edema: number; enhancing: number; total: number };
  status: string;
}

const modalityInfo: Record<ModalityType, { label: string; hint: string }> = {
  T1: { label: 'T1', hint: 'T1-weighted (structural)' },
  T1c: { label: 'T1c', hint: 'T1 + contrast (Gd)' },
  T2: { label: 'T2', hint: 'T2-weighted (fluid bright)' },
  FLAIR: { label: 'FLAIR', hint: 'Fluid-attenuated inversion recovery' },
};

export default function BrainTumorSegmentation() {
  // Health
  const [healthChecking, setHealthChecking] = useState(true);
  const [healthOk, setHealthOk] = useState<boolean | null>(null);

  // Inputs
  const [modalities, setModalities] = useState<Modalities>({ T1: null, T1c: null, T2: null, FLAIR: null });
  const isFormValid = useMemo(() => Object.values(modalities).every(Boolean), [modalities]);

  // Run state
  const [loading, setLoading] = useState(false);
  const [progress, setProgress] = useState(0);
  const [error, setError] = useState<string | null>(null);
  const [lastRunAt, setLastRunAt] = useState<string | null>(null);

  // Results
  const [result, setResult] = useState<SegmentationResult | null>(null);
  const [idx, setIdx] = useState<number>(0);

  // Viewer controls
  const [brightness, setBrightness] = useState(100);
  const [contrast, setContrast] = useState(100);
  const [zoomOpen, setZoomOpen] = useState(false);

  // UI
  const [panel, setPanel] = useState<'inputs' | 'progress' | 'results'>('inputs');

  // Health check
  useEffect(() => {
    (async () => {
      try {
        const r = await axios.get(`${API_BASE_URL}/imaging/brain-tumor/health`);
        setHealthOk(Boolean(r.data?.model_loaded));
      } catch {
        setHealthOk(false);
      } finally {
        setHealthChecking(false);
      }
    })();
  }, []);

  // Simulated smooth progress (UX)
  useEffect(() => {
    if (!loading) return;
    let raf = 0;
    const start = performance.now();
    const duration = 90_000;
    const tick = (t: number) => {
      const x = Math.min((t - start) / duration, 1);
      const eased = 1 - (1 - x) * (1 - x);
      setProgress(Math.min(95, Math.round(eased * 100)));
      if (x < 1) raf = requestAnimationFrame(tick);
    };
    raf = requestAnimationFrame(tick);
    return () => cancelAnimationFrame(raf);
  }, [loading]);

  const handleFile = useCallback((m: ModalityType, f: File | null) => {
    setModalities((prev) => ({ ...prev, [m]: f }));
  }, []);

  const formatBytes = (bytes?: number) => {
    if (!bytes && bytes !== 0) return '';
    const k = 1024; const sizes = ['Bytes','KB','MB','GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return `${parseFloat((bytes / Math.pow(k, i)).toFixed(2))} ${sizes[i]}`;
  };
  const fmtVox = (v: number) => (v < 1000 ? `${v.toFixed(0)}` : `${(v / 1000).toFixed(1)}k`);
  const pct = (p: number, t: number) => (t ? (p / t) * 100 : 0);

  const startSegmentation = useCallback(async () => {
    if (!isFormValid) {
      setError('Please upload T1, T1c, T2 and FLAIR sequences (NIfTI).');
      return;
    }
    setError(null);
    setPanel('progress');
    setLoading(true);
    setResult(null);
    setIdx(0);

    const formData = new FormData();
    (Object.keys(modalities) as ModalityType[]).forEach((m) => {
      const f = modalities[m];
      if (f) formData.append(m, f);
    });

    try {
      const r = await axios.post<SegmentationResult>(`${API_BASE_URL}/imaging/brain-tumor/segment`, formData, {
        headers: { 'Content-Type': 'multipart/form-data' },
        timeout: 300000,
      });
      setResult(r.data);
      setProgress(100);
      setLastRunAt(new Date().toISOString());
      setPanel('results');
      if (r.data.segmented_slices?.length) setIdx(0);
    } catch (err: any) {
      const detail = err?.response?.data?.detail || err?.message || 'Unknown error';
      setError(`Segmentation failed: ${detail}`);
      setPanel('inputs');
    } finally {
      setLoading(false);
    }
  }, [isFormValid, modalities]);

  if (healthChecking) {
    return (
      <div className="min-h-screen flex flex-col items-center justify-center gap-3 p-6">
        <Spinner size="lg" color="primary" />
        <p className="text-gray-600">Connecting to NeuroSeg backend…</p>
      </div>
    );
  }
  if (healthOk === false) {
    return (
      <div className="min-h-screen flex flex-col items-center justify-center gap-4 p-6 text-center">
        <div className="text-6xl">⚠️</div>
        <h1 className="text-2xl font-bold text-red-600">Service Unavailable</h1>
        <p className="max-w-md text-gray-600">The NeuroSeg backend is currently unreachable. Check the service and retry.</p>
        <Button color="primary" onPress={() => location.reload()}>Retry</Button>
      </div>
    );
  }

  const loadedCount = Object.values(modalities).filter(Boolean).length;
  const slices = result?.segmented_slices || [];
  const sliceCount = slices.length;

  return (
    <div className="max-w-7xl mx-auto p-4 md:p-8 bg-gray-50 min-h-screen">
      {/* Sticky KPI Header */}
      <div className="sticky top-0 z-10 mb-6">
        <Card className="border border-gray-200 shadow-sm">
          <CardBody className="py-5">
            <div className="flex flex-col gap-6 md:flex-row md:items-center md:justify-between">
              <div>
                <h1 className="text-2xl md:text-3xl font-semibold tracking-tight text-gray-900">
                  NeuroSeg — MRI Brain Tumor Segmentation
                </h1>
                <p className="text-sm text-gray-600 mt-1">Multimodal glioma segmentation (T1, T1c, T2, FLAIR).</p>
                <div className="mt-2 flex flex-wrap gap-2">
                  <Chip variant="flat" color="primary">v1.8.0</Chip>
                  <Chip variant="flat">3D U-Net + attention</Chip>
                  <Chip variant="flat">NIfTI .nii / .nii.gz</Chip>
                </div>
              </div>
              <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                <Kpi label="Service" value={healthOk ? 'Online' : 'Offline'} tone={healthOk ? 'ok' : 'bad'} />
                <Kpi label="Sequences" value={`${loadedCount}/4`} />
                <Kpi label="Slices" value={sliceCount ? String(sliceCount) : '—'} />
                <Kpi label="Last run" value={lastRunAt ? new Date(lastRunAt).toLocaleString() : '—'} />
              </div>
            </div>
          </CardBody>
        </Card>
      </div>

      {/* Main grid: Left inputs, Center viewer, Right metrics */}
      <div className="grid grid-cols-1 lg:grid-cols-12 gap-6">
        {/* Left: Inputs / Actions */}
        <section className="lg:col-span-3">
          <Card className="border border-gray-200 shadow-sm">
            <CardHeader className="bg-gray-50"><p className="text-lg font-semibold">1 · Input MRI Sequences</p></CardHeader>
            <CardBody>
              <div className="grid grid-cols-1 gap-3">
                {(Object.keys(modalities) as ModalityType[]).map((m) => {
                  const file = modalities[m];
                  return (
                    <Tooltip key={m} content={modalityInfo[m].hint} placement="top">
                      <div className={`rounded-md border ${file ? 'border-emerald-300' : 'border-gray-200'} bg-white p-3`}>
                        <div className="flex items-center justify-between mb-2">
                          <span className="font-medium">{modalityInfo[m].label}</span>
                          {file ? <Chip size="sm" color="success" variant="flat">Loaded</Chip> : <Chip size="sm" variant="flat">Awaiting</Chip>}
                        </div>
                        <input
                          type="file"
                          accept=".nii,.nii.gz"
                          className="block w-full text-sm text-gray-600 file:mr-2 file:py-1 file:px-3 file:rounded-md file:border-0 file:bg-indigo-50 file:text-indigo-700 hover:file:bg-indigo-100"
                          onChange={(e) => handleFile(m, e.target.files?.[0] ?? null)}
                          disabled={loading}
                        />
                        {file && (
                          <div className="mt-1 text-[11px] text-gray-500 break-all">
                            {file.name} • {formatBytes(file.size)}
                          </div>
                        )}
                      </div>
                    </Tooltip>
                  );
                })}
              </div>

              <div className="mt-5 flex flex-col items-stretch gap-2">
                <Button
                  color="primary"
                  isDisabled={!isFormValid || loading}
                  isLoading={loading}
                  onPress={startSegmentation}
                >
                  {loading ? 'Processing…' : 'Run Segmentation'}
                </Button>
                {!isFormValid && <p className="text-xs text-gray-500 text-center">Upload T1, T1c, T2, FLAIR to enable.</p>}
                {error && <div className="rounded-md border-l-4 border-red-600 bg-red-50 p-3 text-sm text-red-800">{error}</div>}
              </div>
            </CardBody>
          </Card>

          {/* Progress */}
          <Card className="mt-6 border border-gray-200 shadow-sm">
            <CardHeader className="bg-gray-50"><p className="text-lg font-semibold">2 · Processing Status</p></CardHeader>
            <CardBody>
              {panel === 'progress' ? (
                <>
                  <div className="flex justify-between text-sm mb-1">
                    <span>Segmentation in progress</span><span>{progress}%</span>
                  </div>
                  <Progress value={progress} color="primary" className="h-2" />
                  <div className="grid grid-cols-2 gap-2 text-xs text-gray-700 mt-4">
                    <Stage label="Pre-processing" done={progress > 25} active={progress > 0 && progress <= 25} />
                    <Stage label="Registration"  done={progress > 50} active={progress > 25 && progress <= 50} />
                    <Stage label="Segmentation"  done={progress > 90} active={progress > 50 && progress <= 90} />
                    <Stage label="Post-processing" done={progress >= 95} active={progress > 90 && progress < 95} />
                  </div>
                </>
              ) : (
                <p className="text-sm text-gray-500">No analysis in progress.</p>
              )}
            </CardBody>
          </Card>
        </section>

        {/* Center: Viewer */}
        <section className="lg:col-span-6">
          <Card className="border border-gray-200 shadow-sm">
            <CardHeader className="bg-gray-50 flex items-center justify-between">
              <p className="font-semibold">3 · Viewer</p>
              <div className="flex items-center gap-3 text-xs text-gray-600">
                <span>WL/WW</span>
                <span className="font-mono">{brightness}%/{contrast}%</span>
                <Button size="sm" variant="flat" onPress={() => setZoomOpen(true)}>Open Lightbox</Button>
              </div>
            </CardHeader>
            <CardBody>
              {!result ? (
                <div className="h-[440px] flex items-center justify-center text-gray-500">Run an analysis to view slices.</div>
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
                        <Slider minValue={50} maxValue={150} value={contrast}  onChange={(v)=>setContrast(Number(v))} aria-label="Contrast"  size="sm" />
                      </div>
                    </div>
                  </div>

                  <Divider className="my-4" />

                  <div className="flex flex-wrap items-center gap-4 text-xs">
                    <Legend swatch="bg-red-500" label="Necrotic / Non-enhancing" />
                    <Legend swatch="bg-blue-500" label="Edema" />
                    <Legend swatch="bg-green-500" label="Enhancing tumor" />
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

        {/* Right: Metrics */}
        <section className="lg:col-span-3">
          <Card className="border border-gray-200 shadow-sm">
            <CardHeader className="bg-gray-50"><p className="font-semibold">4 · Tumor Metrics</p></CardHeader>
            <CardBody>
              {!result ? (
                <p className="text-sm text-gray-500">Metrics will appear after segmentation.</p>
              ) : (
                <>
                  <div className="grid grid-cols-2 gap-3 mb-4">
                    <Metric label="Total (voxels)" value={fmtVox(result.tumor_volumes.total)} tone="base" />
                    <Metric label="Necrotic" value={fmtVox(result.tumor_volumes.necrotic)} tone="red" />
                    <Metric label="Edema" value={fmtVox(result.tumor_volumes.edema)} tone="blue" />
                    <Metric label="Enhancing" value={fmtVox(result.tumor_volumes.enhancing)} tone="green" />
                  </div>

                  {result.tumor_volumes.total > 0 && (
                    <>
                      <p className="text-xs text-gray-600 mb-1">Composition</p>
                      <div className="w-full h-3 rounded-full overflow-hidden flex">
                        <div className="bg-red-500"   style={{ width: `${pct(result.tumor_volumes.necrotic,  result.tumor_volumes.total)}%` }} />
                        <div className="bg-blue-500"  style={{ width: `${pct(result.tumor_volumes.edema,     result.tumor_volumes.total)}%` }} />
                        <div className="bg-green-500" style={{ width: `${pct(result.tumor_volumes.enhancing, result.tumor_volumes.total)}%` }} />
                      </div>
                      <div className="flex justify-between text-[11px] text-gray-500 mt-1">
                        <span>{Math.round(pct(result.tumor_volumes.necrotic,  result.tumor_volumes.total))}% Necrotic</span>
                        <span>{Math.round(pct(result.tumor_volumes.edema,     result.tumor_volumes.total))}% Edema</span>
                        <span>{Math.round(pct(result.tumor_volumes.enhancing, result.tumor_volumes.total))}% Enhancing</span>
                      </div>
                    </>
                  )}
                </>
              )}
            </CardBody>
            {result?.download_url && (
              <CardFooter className="bg-gray-50">
                <Button as="a" href={result.download_url} color="primary" variant="bordered" className="w-full">
                  Download Full Segmentation
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
                alt="Zoomed segmented slice"
                className="w-full h-full object-contain"
                style={{ filter: `brightness(${brightness}%) contrast(${contrast}%)` }}
              />
            )}
          </ModalBody>
          <ModalFooter className="justify-between">
            <div className="flex-1 pr-6">
              <div className="grid grid-cols-2 gap-3">
                <div>
                  <p className="text-xs text-gray-600 mb-1">Brightness</p>
                  <Slider minValue={50} maxValue={150} value={brightness} onChange={(v)=>setBrightness(Number(v))} size="sm" />
                </div>
                <div>
                  <p className="text-xs text-gray-600 mb-1">Contrast</p>
                  <Slider minValue={50} maxValue={150} value={contrast} onChange={(v)=>setContrast(Number(v))} size="sm" />
                </div>
              </div>
            </div>
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
function Stage({ label, done, active }: { label: string; done?: boolean; active?: boolean }) {
  return (
    <div className="rounded bg-gray-50 p-2">
      <div className="flex items-center gap-2">
        <span className={`inline-block h-2 w-2 rounded-full ${done ? 'bg-emerald-500' : active ? 'bg-amber-500' : 'bg-gray-300'}`} />
        <span className="text-[12px]">{label}</span>
      </div>
      <div className="ml-4 text-[11px] text-gray-500">
        {done ? 'Complete' : active ? 'In progress' : 'Pending'}
      </div>
    </div>
  );
}
function Metric({ label, value, tone = 'base' }: { label: string; value: string; tone?: 'base'|'red'|'blue'|'green' }) {
  const map: Record<string, string> = {
    base: 'text-gray-900',
    red: 'text-red-600',
    blue: 'text-blue-600',
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
