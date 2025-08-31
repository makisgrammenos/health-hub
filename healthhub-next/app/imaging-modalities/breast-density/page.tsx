'use client';

import React, { useRef, useState, useMemo } from 'react';
import axios from 'axios';
import { Card, CardBody, CardFooter, CardHeader } from '@nextui-org/card';
import { Button } from '@nextui-org/button';
import { Progress } from '@nextui-org/progress';
import { Chip } from '@nextui-org/chip';
import { Tooltip } from '@nextui-org/tooltip';
import { Table, TableHeader, TableBody, TableColumn, TableRow, TableCell } from '@nextui-org/table';
import { Modal, ModalBody, ModalContent, ModalFooter } from '@nextui-org/modal';
import { Slider } from '@nextui-org/slider';
import { Divider } from '@nextui-org/divider';
import { Spinner } from '@nextui-org/spinner';
import { Download, Upload, Check, X, Trash2, Info } from 'lucide-react';

const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

type ProbMap = Record<string, number>;

interface PredictionResponse {
  filename: string;
  prediction_index: number;
  prediction_label: string;   // e.g. A/B/C/D or I/II/III/IV
  probabilities: ProbMap;
}

const ACCEPTED_MIME = new Set([
  'image/jpeg',
  'image/png',
  'image/tiff',
  'image/bmp',
]);

const LABEL_META: Record<string, { short: string; long: string; tone: 'success'|'primary'|'warning'|'danger' }> = {
  A: { short: 'A • Fatty', long: 'Almost entirely fatty (low fibroglandular tissue).', tone: 'success' },
  B: { short: 'B • Scattered FG', long: 'Scattered fibroglandular densities (mild).', tone: 'primary' },
  C: { short: 'C • Heterogeneous', long: 'Heterogeneously dense; may obscure small masses.', tone: 'warning' },
  D: { short: 'D • Extremely dense', long: 'Extremely dense; reduces mammographic sensitivity.', tone: 'danger' },
  I: { short: 'I • Fatty', long: 'Almost entirely fatty (alternate scale).', tone: 'success' },
  II: { short: 'II • Scattered FG', long: 'Scattered fibroglandular densities.', tone: 'primary' },
  III: { short: 'III • Heterogeneous', long: 'Heterogeneously dense.', tone: 'warning' },
  IV: { short: 'IV • Extremely dense', long: 'Extremely dense; lowers sensitivity.', tone: 'danger' },
};

function toneFor(label: string): 'success'|'primary'|'warning'|'danger' {
  const L = label.trim().toUpperCase();
  if (LABEL_META[L]) return LABEL_META[L].tone;
  if (['A','I'].includes(L)) return 'success';
  if (['B','II'].includes(L)) return 'primary';
  if (['C','III'].includes(L)) return 'warning';
  return 'danger';
}

export default function BreastDensityProfessional() {
  const fileInputRef = useRef<HTMLInputElement>(null);

  // State
  const [file, setFile] = useState<File | null>(null);
  const [preview, setPreview] = useState<string | null>(null);

  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<PredictionResponse | null>(null);

  // Viewer controls
  const [zoomOpen, setZoomOpen] = useState(false);
  const [brightness, setBrightness] = useState(100);
  const [contrast, setContrast] = useState(100);

  const isReadyToRun = !!file && !loading;

  const sortedProbs = useMemo(() => {
    if (!result?.probabilities) return [];
    return Object.entries(result.probabilities).sort((a, b) => b[1] - a[1]);
  }, [result]);

  // Handlers
  const openPicker = () => fileInputRef.current?.click();

  const onFilePicked = (f: File) => {
    if (!ACCEPTED_MIME.has(f.type)) {
      setError(`Unsupported file type: ${f.type || 'unknown'}. Use JPEG/PNG/TIFF/BMP.`);
      setFile(null);
      setPreview(null);
      setResult(null);
      return;
    }
    setError(null);
    setResult(null);
    setFile(f);
    setPreview(URL.createObjectURL(f));
  };

  const onInputChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const f = e.target.files?.[0];
    if (f) onFilePicked(f);
  };

  const clearSelection = () => {
    setFile(null);
    setPreview(null);
    setResult(null);
    setError(null);
    setBrightness(100);
    setContrast(100);
    if (fileInputRef.current) fileInputRef.current.value = '';
  };

  const onDrop: React.DragEventHandler<HTMLDivElement> = (e) => {
    e.preventDefault();
    const f = e.dataTransfer.files?.[0];
    if (f) onFilePicked(f);
  };

  const runInference = async () => {
    if (!file) return;
    setLoading(true);
    setError(null);
    setResult(null);
    const form = new FormData();
    form.append('file', file);

    try {
      const { data } = await axios.post<PredictionResponse>(`${API_BASE_URL}/imaging/breast-density/predict`, form, {
        headers: { 'Content-Type': 'multipart/form-data' },
        timeout: 120000,
      });
      console.log(data)
      setResult(data);
    } catch (err: any) {
      const detail = err?.response?.data?.detail || err?.message || 'Unknown error';
      setError(`Prediction failed: ${detail}`);
    } finally {
      setLoading(false);
    }
  };

  const exportJSON = () => {
    if (!result) return;
    const name = (result.filename || 'breast_density').replace(/\s+/g, '_');
    const blob = new Blob([JSON.stringify(result, null, 2)], { type: 'application/json;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = `${name}.json`; a.click();
    URL.revokeObjectURL(url);
  };

  return (
    <div className="mx-auto max-w-7xl p-4 md:p-8 bg-gray-50 min-h-screen">
      {/* Header / KPIs */}
      <Card className="mb-8 border border-gray-200 shadow-sm">
        <CardBody className="py-5">
          <div className="flex flex-col gap-6 md:flex-row md:items-center md:justify-between">
            <div>
              <h1 className="text-3xl font-semibold tracking-tight text-gray-900">
                BI-RADS Breast Density – Professional
              </h1>
              <p className="text-sm text-gray-600 mt-1">
                Single-image classifier returning BI-RADS density with calibrated per-class probabilities.
              </p>
              <div className="mt-2 flex flex-wrap gap-2">
                <Chip variant="flat" color="primary">Model: BiradsDensityModel</Chip>
                <Chip variant="flat">Inputs: JPEG · PNG · TIFF · BMP</Chip>
                <Chip variant="flat">Endpoint: POST /predict</Chip>
              </div>
            </div>

            <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
              <Kpi label="Status" value="Online" tone="ok" />
              <Kpi label="File" value={file ? 'Selected' : '—'} />
              <Kpi label="Result" value={result ? 'Available' : '—'} />
              <Kpi label="Top Confidence" value={result ? `${Math.round((sortedProbs[0]?.[1] || 0) * 100)}%` : '—'} />
            </div>
          </div>
        </CardBody>
      </Card>

      {/* Main layout */}
      <div className="grid grid-cols-1 lg:grid-cols-12 gap-6">
        {/* Left: Upload & Actions */}
        <section className="lg:col-span-4">
          <Card className="border border-gray-200 shadow-sm">
            <CardHeader className="bg-gray-50">
              <p className="text-lg font-semibold">1 · Input Mammogram</p>
            </CardHeader>
            <CardBody>
              <div
                onDrop={onDrop}
                onDragOver={(e) => e.preventDefault()}
                className="
                  border-2 border-dashed rounded-lg p-6 bg-white
                  border-gray-300 hover:border-gray-400 transition
                  text-center
                "
                aria-label="Upload area"
              >
                <input
                  ref={fileInputRef}
                  type="file"
                  accept="image/jpeg,image/png,image/tiff,image/bmp"
                  onChange={onInputChange}
                  className="hidden"
                />
                <Upload className="mx-auto mb-2 text-gray-400" size={40} />
                <p className="text-gray-800">Drag & drop an image here</p>
                <p className="text-sm text-gray-500 mb-4">or</p>
                <div className="flex items-center justify-center gap-2">
                  <Button color="primary" variant="flat" onPress={openPicker}>Browse Files</Button>
                  {file && (
                    <Button color="danger" variant="light" onPress={clearSelection} startContent={<Trash2 size={16} />}>
                      Clear
                    </Button>
                  )}
                </div>
                {file && (
                  <div className="mt-3 text-sm text-emerald-700 flex items-center justify-center">
                    <Check size={16} className="mr-1" /> {file.name}
                  </div>
                )}
              </div>

              {error && (
                <div className="mt-4 rounded-md border-l-4 border-red-600 bg-red-50 p-3 text-sm text-red-800 flex items-start gap-2">
                  <X size={18} className="mt-0.5 shrink-0" />
                  <span>{error}</span>
                </div>
              )}

              <div className="mt-6 flex flex-col gap-2">
                <Button
                  color="primary"
                  size="lg"
                  isDisabled={!isReadyToRun}
                  isLoading={loading}
                  onPress={runInference}
                  className="w-full"
                >
                  {loading ? 'Analyzing…' : 'Run Classification'}
                </Button>
                {!file && (
                  <p className="text-xs text-gray-500 text-center">
                    Tip: JPEG/PNG/TIFF/BMP up to typical mammogram sizes. PHI must be removed before upload.
                  </p>
                )}
              </div>
            </CardBody>
          </Card>

          {/* Viewer controls (sticky on smaller screens) */}
          <Card className="mt-6 border border-gray-200 shadow-sm">
            <CardHeader className="bg-gray-50 flex items-center justify-between">
              <p className="font-semibold">Viewer Controls</p>
              {preview && (
                <Button size="sm" variant="flat" onPress={() => setZoomOpen(true)}>
                  Open Lightbox
                </Button>
              )}
            </CardHeader>
            <CardBody>
              {!preview ? (
                <div className="text-sm text-gray-500">Load an image to adjust brightness/contrast.</div>
              ) : (
                <div className="grid grid-cols-2 gap-4">
                  <div>
                    <div className="flex items-center justify-between text-xs text-gray-600 mb-1">
                      <span>Brightness</span><span>{brightness}%</span>
                    </div>
                    <Slider minValue={50} maxValue={150} value={brightness} onChange={(v)=>setBrightness(Number(v))} size="sm" />
                  </div>
                  <div>
                    <div className="flex items-center justify-between text-xs text-gray-600 mb-1">
                      <span>Contrast</span><span>{contrast}%</span>
                    </div>
                    <Slider minValue={50} maxValue={150} value={contrast} onChange={(v)=>setContrast(Number(v))} size="sm" />
                  </div>
                </div>
              )}
            </CardBody>
          </Card>
        </section>

        {/* Center: Image Viewer */}
        <section className="lg:col-span-4">
          <Card className="border border-gray-200 shadow-sm">
            <CardHeader className="bg-gray-50">
              <p className="font-semibold">2 · Image Viewer</p>
            </CardHeader>
            <CardBody>
              {!preview ? (
                <div className="h-[460px] flex items-center justify-center text-gray-500">
                  No image selected.
                </div>
              ) : (
                <div className="relative w-full h-[460px] rounded-md overflow-hidden bg-black">
                  {/* eslint-disable-next-line @next/next/no-img-element */}
                  <img
                    src={preview}
                    alt="Uploaded mammogram"
                    className="h-full w-full object-contain"
                    style={{ filter: `brightness(${brightness}%) contrast(${contrast}%)` }}
                  />
                </div>
              )}
            </CardBody>
          </Card>
        </section>

        {/* Right: Results */}
        <section className="lg:col-span-4">
          <Card className="border border-gray-200 shadow-sm">
            <CardHeader className="bg-gray-50 flex items-center justify-between">
              <p className="font-semibold">3 · Classification Result</p>
              <Tooltip content="BI-RADS density is a qualitative category; use with clinical context.">
                <Info size={16} className="text-gray-500" />
              </Tooltip>
            </CardHeader>
            <CardBody>
              {loading && (
                <div className="h-[200px] flex flex-col items-center justify-center gap-3">
                  <Spinner color="primary" />
                  <p className="text-sm text-gray-600">Running inference…</p>
                </div>
              )}

              {!loading && !result && (
                <div className="h-[200px] flex items-center justify-center text-gray-500">
                  Results will appear here after analysis.
                </div>
              )}

              {!loading && result && (
                <>
                  <div className="flex flex-col items-center text-center mb-4">
                    <Chip size="lg" color={toneFor(result.prediction_label)} className="mb-2">
                      {LABEL_META[result.prediction_label]?.short || `Density ${result.prediction_label}`}
                    </Chip>
                    <p className="text-sm text-gray-600">
                      {LABEL_META[result.prediction_label]?.long || 'Breast density category per BI-RADS.'}
                    </p>
                    <p className="text-xs text-gray-500">File: {result.filename}</p>
                  </div>

                  <Divider className="my-4" />

                  <p className="text-sm font-medium mb-2">Per-class probabilities</p>
                  <Table aria-label="Per-class probabilities">
                    <TableHeader>
                      <TableColumn width={120}>CLASS</TableColumn>
                      <TableColumn>CONFIDENCE</TableColumn>
                      <TableColumn width={80} className="text-right">VALUE</TableColumn>
                    </TableHeader>
                    <TableBody>
                      {sortedProbs.map(([label, p]) => (
                        <TableRow key={label}>
                          <TableCell>
                            <Chip size="sm" variant="flat" color={toneFor(label)}>{label}</Chip>
                          </TableCell>
                          <TableCell>
                            <Progress value={p * 100} color={toneFor(label)} className="max-w-full" />
                          </TableCell>
                          <TableCell className="text-right text-sm">{(p * 100).toFixed(1)}%</TableCell>
                        </TableRow>
                      ))}
                    </TableBody>
                  </Table>

                  <div className="mt-6 rounded-lg border-l-4 border-blue-600 bg-gray-50 p-4">
                    <p className="text-sm font-semibold mb-1">Clinical Note</p>
                    <p className="text-sm text-gray-700">
                      Breast density influences mammographic sensitivity. Combine with radiologist assessment and patient history.
                    </p>
                  </div>
                </>
              )}
            </CardBody>
            <CardFooter className="bg-gray-50">
              <div className="w-full flex items-center justify-between">
                <p className="text-xs text-gray-500">
                  Decision support only; not a standalone diagnostic.
                </p>
                <Button
                  size="sm"
                  variant="flat"
                  color="primary"
                  isDisabled={!result}
                  onPress={exportJSON}
                  endContent={<Download size={16} />}
                >
                  Export JSON
                </Button>
              </div>
            </CardFooter>
          </Card>
        </section>
      </div>

      {/* Lightbox */}
      <Modal isOpen={zoomOpen} onClose={() => setZoomOpen(false)} size="5xl">
        <ModalContent>
          <ModalBody className="p-0 bg-black">
            {!preview ? (
              <div className="h-[70vh] w-full flex items-center justify-center text-gray-400">
                No image selected.
              </div>
            ) : (
              // eslint-disable-next-line @next/next/no-img-element
              <img
                src={preview}
                alt="Preview zoom"
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

/* ---------- UI helpers ---------- */
function Kpi({ label, value, tone }: { label: string; value: string; tone?: 'ok'|'bad'|'base' }) {
  const color = tone === 'ok' ? 'text-emerald-700' : tone === 'bad' ? 'text-red-700' : 'text-gray-900';
  return (
    <div className="rounded-lg border border-gray-200 bg-white p-3 text-center">
      <p className="text-[11px] uppercase tracking-wide text-gray-500">{label}</p>
      <p className={`font-mono text-xl font-semibold ${color}`}>{value}</p>
    </div>
  );
}
