'use client';

import React, { useState, useCallback, useEffect } from 'react';
import axios from 'axios';
import { Card, CardBody, CardFooter } from '@nextui-org/card';
import { Button } from '@nextui-org/button';
import { Progress } from '@nextui-org/progress';
import { Tooltip } from '@nextui-org/tooltip';
import { Chip } from '@nextui-org/chip';
import { Image } from '@nextui-org/image';
import { Modal, ModalContent, ModalBody, ModalFooter } from '@nextui-org/modal';
import { Skeleton } from '@nextui-org/skeleton';
import { Spinner } from '@nextui-org/spinner';
import  {Slider} from '@nextui-org/slider';

// Define the API base URL for production/development environments
const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

interface SegmentationSlice {
  slice_index: number;
  image: string;
}

interface SegmentationResult {
  segmented_slices: SegmentationSlice[];
  download_url: string;
  tumor_volumes: {
    necrotic: number;
    edema: number;
    enhancing: number;
    total: number;
  };
  status: string;
}

type ModalityType = 'T1c' | 'T1' | 'T2' | 'FLAIR';

type Modalities = {
  [key in ModalityType]: File | null;
};

const modalityDescriptions = {
  'T1c': 'T1-weighted with contrast enhancement (T1c)',
  'T1': 'T1-weighted MRI',
  'T2': 'T2-weighted MRI',
  'FLAIR': 'Fluid Attenuated Inversion Recovery'
};

export default function BrainTumorSegmentation() {
  const [modalities, setModalities] = useState<Modalities>({
    T1c: null,
    T1: null,
    T2: null,
    FLAIR: null,
  });
  const [segmentationResult, setSegmentationResult] = useState<SegmentationResult | null>(null);
  const [loading, setLoading] = useState(false);
  const [healthStatus, setHealthStatus] = useState<boolean | null>(null);
  const [healthChecking, setHealthChecking] = useState(true);
  const [progress, setProgress] = useState(0);
  const [error, setError] = useState<string | null>(null);
  const [selectedSliceIndex, setSelectedSliceIndex] = useState<number | null>(null);
  const [imageModalOpen, setImageModalOpen] = useState(false);
  const [currentImageUrl, setCurrentImageUrl] = useState<string | null>(null);

  // Check API health on component mount
  useEffect(() => {
    const checkHealth = async () => {
      try {
        const response = await axios.get(`http://localhost:8000/imaging/brain-tumor/health`);
        setHealthStatus(response.data.model_loaded);
      } catch (err) {
        setHealthStatus(false);
      } finally {
        setHealthChecking(false);
      }
    };

    checkHealth();
  }, []);

  const handleFileChange = useCallback((modality: ModalityType, file: File | null) => {
    setModalities((prev) => ({ ...prev, [modality]: file }));
  }, []);

  const isFormValid = Object.values(modalities).every(file => file !== null);

  const simulateProgress = useCallback(() => {
  setProgress(0);
  const startTime = Date.now();
  const simDuration = 120 * 1000; // ~1 minute
  let rafId: number;

  const tick = () => {
    const elapsed = Date.now() - startTime;
    const t = Math.min(elapsed / simDuration, 1);          // 0 → 1 over simDuration
    const eased = 1 - (1 - t) * (1 - t);                    // ease-out quad
    const next = Math.min(eased * 95, 95);                  // cap at 95%
    setProgress(next);

    if (elapsed < simDuration) {
      rafId = requestAnimationFrame(tick);
    }
  };

  rafId = requestAnimationFrame(tick);
  return () => cancelAnimationFrame(rafId);
}, []);


  const startSegmentation = useCallback(async () => {
    if (!isFormValid) {
      setError('Please upload all four MRI modalities to proceed with segmentation analysis.');
      return;
    }

    setLoading(true);
    setError(null);
    const cleanupProgressSimulation = simulateProgress();

    const formData = new FormData();
    Object.entries(modalities).forEach(([key, file]) => {
      if (file) formData.append(key, file);
    });

    try {
      const response = await axios.post(`${API_BASE_URL}/imaging/brain-tumor/segment`, formData, {
        headers: { 'Content-Type': 'multipart/form-data' },
        timeout: 300000, // 5-minute timeout for large files
      });

      setSegmentationResult(response.data);
      setProgress(100);

      // Set the first slice as selected by default if we have results
      if (response.data.segmented_slices && response.data.segmented_slices.length > 0) {
        setSelectedSliceIndex(0);
      }
    } catch (err) {
      const errorMessage = axios.isAxiosError(err) && err.response?.data?.detail
        ? err.response.data.detail
        : 'An error occurred during segmentation analysis. Please verify your input files and try again.';

      setError(errorMessage);
    } finally {
      cleanupProgressSimulation();
      setLoading(false);
    }
  }, [modalities, isFormValid, simulateProgress]);

  const formatBytes = (bytes: number, decimals = 2) => {
    if (bytes === 0) return '0 Bytes';
    const k = 1024;
    const dm = decimals < 0 ? 0 : decimals;
    const sizes = ['Bytes', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return parseFloat((bytes / Math.pow(k, i)).toFixed(dm)) + ' ' + sizes[i];
  };

  const openImageModal = (imageUrl: string) => {
    setCurrentImageUrl(imageUrl);
    setImageModalOpen(true);
  };

  // Function to format voxel volumes
  const formatVoxelVolume = (voxels: number) => {
    if (voxels === 0) return '0';
    if (voxels < 1000) return voxels.toFixed(0);
    return `${(voxels / 1000).toFixed(1)}k`;
  };

  // Calculate percentage of each tumor component
  const calculatePercentage = (component: number, total: number) => {
    if (total === 0) return 0;
    return (component / total) * 100;
  };

  if (healthChecking) {
    return (
      <div className="flex flex-col items-center justify-center min-h-screen p-4">
        <Spinner size="lg" color="primary" />
        <p className="mt-4 text-gray-600">Connecting to NeuroSeg backend...</p>
      </div>
    );
  }

  if (healthStatus === false) {
    return (
      <div className="flex flex-col items-center justify-center min-h-screen p-4 text-center">
        <div className="text-red-600 text-6xl mb-4">⚠️</div>
        <h1 className="text-2xl font-bold text-red-600 mb-2">Service Unavailable</h1>
        <p className="mb-6 text-gray-600">
          The NeuroSeg backend service is currently unavailable. Please try again later or contact support.
        </p>
        <Button
          color="primary"
          onClick={() => window.location.reload()}
        >
          Retry Connection
        </Button>
      </div>
    );
  }

  return (
    <div className="max-w-7xl mx-auto p-4 md:p-8 bg-gray-50 min-h-screen">
      <div className="mb-8 text-center">
        <h1 className="text-3xl md:text-4xl font-bold mb-4">
          <span className="text-blue-700">NeuroSeg</span>: Brain Tumor Segmentation Analysis
        </h1>
        <div className="max-w-3xl mx-auto">
          <p className="text-gray-700">
            Advanced multimodal MRI-based glioma segmentation using convolutional neural networks.
            Upload standard MRI sequences to automatically identify and segment tumor regions.
          </p>
        </div>
      </div>

      {/* Step 1: Upload Files */}
      <Card className="mb-8 shadow-md">
        <CardBody>
          <h2 className="text-xl font-bold mb-4 text-blue-800 flex items-center">
            <span className="bg-blue-100 text-blue-800 rounded-full w-8 h-8 inline-flex items-center justify-center mr-2">1</span>
            Input MRI Sequences
          </h2>
          <p className="mb-6 text-gray-600 text-sm">
            Upload the standard brain MRI sequences (NIfTI format: .nii or .nii.gz). All four modalities are required for accurate segmentation.
          </p>
          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-4">
            {Object.keys(modalities).map((modality) => {
              const mod = modality as ModalityType;
              const file = modalities[mod];
              return (
                <Tooltip key={mod} content={modalityDescriptions[mod]}>
                  <Card shadow="sm" className={`border ${file ? 'border-green-300' : 'border-gray-200'}`}>
                    <CardBody className="p-4">
                      <div className="flex justify-between items-center mb-2">
                        <h3 className="font-semibold">{mod}</h3>
                        {file && (
                          <Chip size="sm" color="success" variant="flat">
                            Loaded
                          </Chip>
                        )}
                      </div>
                      <input
                        type="file"
                        accept=".nii,.nii.gz"
                        onChange={(e) => handleFileChange(mod, e.target.files?.[0] || null)}
                        className="block w-full text-sm text-gray-500 file:mr-2 file:py-1 file:px-3 file:rounded-md file:border-0 file:bg-blue-50 file:text-blue-700 hover:file:bg-blue-100"
                        disabled={loading}
                      />
                      {file && (
                        <div className="mt-2 text-xs text-gray-500">
                          {file.name.length > 20 ? `${file.name.substring(0, 20)}...` : file.name} ({formatBytes(file.size)})
                        </div>
                      )}
                    </CardBody>
                  </Card>
                </Tooltip>
              );
            })}
          </div>
          <div className="text-center mt-6">
            <Button
              isDisabled={!isFormValid || loading}
              isLoading={loading}
              onPress={startSegmentation}
              color="primary"
              className="bg-blue-700 text-white font-medium px-6"
              size="lg"
              // startContent={!loading && <span className="material-icons text-sm">science</span>}
            >
              {loading ? "Processing..." : "Execute Segmentation Analysis"}
            </Button>
          </div>
          {error && (
            <div className="mt-4 p-3 bg-red-50 border border-red-200 text-red-700 rounded-md text-center">
              {error}
            </div>
          )}
        </CardBody>
      </Card>

      {/* Step 2: Segmentation Progress */}
      <Card className="mb-8 shadow-md">
        <CardBody>
          <h2 className="text-xl font-bold mb-4 text-blue-800 flex items-center">
            <span className="bg-blue-100 text-blue-800 rounded-full w-8 h-8 inline-flex items-center justify-center mr-2">2</span>
            Processing Status
          </h2>
          {loading ? (
            <div>
              <div className="flex justify-between mb-2">
                <span className="text-sm font-medium text-gray-700">Segmentation in progress</span>
                <span className="text-sm font-medium text-gray-700">{Math.round(progress)}%</span>
              </div>
              <Progress
                color="primary"
                value={progress}
                className="h-2"
                aria-label="Segmentation progress"
              />
              <div className="mt-4 grid grid-cols-2 md:grid-cols-4 gap-2 text-xs text-gray-600">
                <div className="bg-gray-50 p-2 rounded">
                  <div>Pre-processing</div>
                  <div className={progress > 25 ? "text-green-600 font-medium" : ""}>
                    {progress > 25 ? "Complete" : "Pending"}
                  </div>
                </div>
                <div className="bg-gray-50 p-2 rounded">
                  <div>Image Registration</div>
                  <div className={progress > 50 ? "text-green-600 font-medium" : ""}>
                    {progress > 30 ? "Complete" : progress > 25 ? "In Progress" : "Pending"}
                  </div>
                </div>
                <div className="bg-gray-50 p-2 rounded">
                  <div>Segmentation</div>
                  <div className={progress > 75 ? "text-green-600 font-medium" : ""}>
                    {progress > 90 ? "Complete" : progress > 30 ? "In Progress" : "Pending"}
                  </div>
                </div>
                <div className="bg-gray-50 p-2 rounded">
                  <div>Post-processing</div>
                  <div className={progress > 95 ? "text-green-600 font-medium" : ""}>
                    {progress > 95 ? "Complete" : progress > 90 ? "In Progress" : "Pending"}
                  </div>
                </div>
              </div>
            </div>
          ) : (
            <div className="p-6 text-center text-gray-500">
              No analysis in progress. Upload MRI sequences and press "Execute Segmentation Analysis" to begin.
            </div>
          )}
        </CardBody>
      </Card>

      {/* Step 3: Segmentation Results */}
      <Card className="shadow-md">
        <CardBody>
          <h2 className="text-xl font-bold mb-4 text-blue-800 flex items-center">
            <span className="bg-blue-100 text-blue-800 rounded-full w-8 h-8 inline-flex items-center justify-center mr-2">3</span>
            Segmentation Results
          </h2>
          {segmentationResult ? (
            <div>
              <div className="mb-6 bg-blue-50 p-4 rounded-lg">
                <h3 className="font-medium text-blue-800 mb-2">Analysis Complete</h3>
                <p className="text-sm text-gray-700">
                  The segmentation identified tumor regions across {segmentationResult.segmented_slices.length} key MRI slices.
                  Each image shows the original MRI with color overlay for different tumor components.
                  View results below or download the full volumetric segmentation for further analysis.
                </p>
              </div>

              <div className="mb-6">
                <h3 className="font-medium mb-2">Slice Selection</h3>
                <div className="flex overflow-x-auto pb-2 space-x-2">
                  {segmentationResult.segmented_slices.map((slice, idx) => (
                    <Button
                      key={idx}
                      size="sm"
                      variant={selectedSliceIndex === idx ? "solid" : "flat"}
                      color={selectedSliceIndex === idx ? "primary" : "default"}
                      onPress={() => setSelectedSliceIndex(idx)}
                    >
                      Slice {slice.slice_index}
                    </Button>
                  ))}
                </div>
              </div>

              {selectedSliceIndex !== null && segmentationResult.segmented_slices[selectedSliceIndex] ? (
                <div className="flex flex-col items-center mb-6">
                  <Card className="w-full max-w-2xl shadow-md">
                    <CardBody className="p-0 relative">
                      <Image
                        src={segmentationResult.segmented_slices[selectedSliceIndex].image}
                        alt={`Segmented MRI Slice ${segmentationResult.segmented_slices[selectedSliceIndex].slice_index}`}
                        className="w-full h-auto rounded-lg cursor-pointer"
                        radius="lg"
                        onClick={() => openImageModal(segmentationResult.segmented_slices[selectedSliceIndex].image)}
                      />
                      <div className="absolute bottom-2 right-2">
                        <Button
                          size="sm"
                          isIconOnly
                          color="default"
                          variant="flat"
                          className="bg-white/70"
                          onPress={() => openImageModal(segmentationResult.segmented_slices[selectedSliceIndex].image)}
                        >
                          <span className="material-icons text-sm">zoom_in</span>
                        </Button>
                      </div>
                    </CardBody>
                    <CardFooter className="flex flex-col sm:flex-row justify-between items-center gap-2">
                      <div className="text-sm">
                        <span className="font-medium">Slice {segmentationResult.segmented_slices[selectedSliceIndex].slice_index}</span>
                        <div className="text-xs text-gray-500">Axial View with Segmentation Overlay</div>
                      </div>
                      <div className="flex flex-wrap items-center gap-2">
                        <div className="flex items-center">
                          <span className="w-3 h-3 bg-red-500 rounded-full mr-1"></span>
                          <span className="text-xs">Necrotic/Non-enhancing</span>
                        </div>
                        <div className="flex items-center">
                          <span className="w-3 h-3 bg-blue-500 rounded-full mr-1"></span>
                          <span className="text-xs">Edema</span>
                        </div>
                        <div className="flex items-center">
                          <span className="w-3 h-3 bg-green-500 rounded-full mr-1"></span>
                          <span className="text-xs">Enhancing Tumor</span>
                        </div>
                      </div>
                    </CardFooter>
                  </Card>
                </div>
              ) : (
                <div className="text-center p-6 bg-gray-50 rounded-lg mb-6">
                  <p className="text-gray-600">Select a slice number above to view detailed segmentation results.</p>
                </div>
              )}

              {/* Tumor Composition Analysis */}
              {segmentationResult.tumor_volumes && (
                <div className="mb-6">
                  <h3 className="font-medium mb-4">Tumor Composition Analysis</h3>
                  <div className="bg-white p-4 rounded-lg shadow-sm">
                    <div className="grid grid-cols-1 md:grid-cols-4 gap-4 mb-4">
                      <div className="bg-gray-50 p-4 rounded-lg text-center">
                        <div className="text-2xl font-bold text-blue-800">
                          {formatVoxelVolume(segmentationResult.tumor_volumes.total)}
                        </div>
                        <div className="text-xs text-gray-600">Total Tumor Volume (voxels)</div>
                      </div>
                      <div className="bg-gray-50 p-4 rounded-lg text-center">
                        <div className="text-2xl font-bold text-red-600">
                          {formatVoxelVolume(segmentationResult.tumor_volumes.necrotic)}
                        </div>
                        <div className="text-xs text-gray-600">Necrotic/Non-enhancing</div>
                      </div>
                      <div className="bg-gray-50 p-4 rounded-lg text-center">
                        <div className="text-2xl font-bold text-blue-600">
                          {formatVoxelVolume(segmentationResult.tumor_volumes.edema)}
                        </div>
                        <div className="text-xs text-gray-600">Edema</div>
                      </div>
                      <div className="bg-gray-50 p-4 rounded-lg text-center">
                        <div className="text-2xl font-bold text-green-600">
                          {formatVoxelVolume(segmentationResult.tumor_volumes.enhancing)}
                        </div>
                        <div className="text-xs text-gray-600">Enhancing Tumor</div>
                      </div>
                    </div>

                    {/* Composition Bar */}
                    {segmentationResult.tumor_volumes.total > 0 && (
                      <div className="mt-4">
                        <div className="text-sm font-medium mb-1">Tumor Composition</div>
                        <div className="w-full h-4 flex rounded-full overflow-hidden">
                          <div
                            className="bg-red-500"
                            style={{ width: `${calculatePercentage(segmentationResult.tumor_volumes.necrotic, segmentationResult.tumor_volumes.total)}%` }}
                          ></div>
                          <div
                            className="bg-blue-500"
                            style={{ width: `${calculatePercentage(segmentationResult.tumor_volumes.edema, segmentationResult.tumor_volumes.total)}%` }}
                          ></div>
                          <div
                            className="bg-green-500"
                            style={{ width: `${calculatePercentage(segmentationResult.tumor_volumes.enhancing, segmentationResult.tumor_volumes.total)}%` }}
                          ></div>
                        </div>
                        <div className="flex justify-between text-xs text-gray-500 mt-1">
                          <span>{Math.round(calculatePercentage(segmentationResult.tumor_volumes.necrotic, segmentationResult.tumor_volumes.total))}% Necrotic</span>
                          <span>{Math.round(calculatePercentage(segmentationResult.tumor_volumes.edema, segmentationResult.tumor_volumes.total))}% Edema</span>
                          <span>{Math.round(calculatePercentage(segmentationResult.tumor_volumes.enhancing, segmentationResult.tumor_volumes.total))}% Enhancing</span>
                        </div>
                      </div>
                    )}
                  </div>
                </div>
              )}
            )}

              {/* Download Button */}
              <div className="text-center mt-4">
                <Button
                  as="a"
                  href={segmentationResult.download_url}
                  color="primary"
                  variant="bordered"
                  startContent={<span className="material-icons text-sm">download</span>}
                >
                  Download Full Segmentation Results
                </Button>
              </div>
            </div>
          ) : (
            <div className="p-6 text-center text-gray-500">
              No segmentation results yet. Complete the analysis to view detailed findings.
            </div>
          )}
        </CardBody>
      </Card>

      {/* Image Modal */}
      <Modal isOpen={imageModalOpen} onClose={() => setImageModalOpen(false)} size="5xl">
        <ModalContent>
          <ModalBody className="p-0">
            {currentImageUrl && (
              <Image
                src={currentImageUrl}
                alt="Full-size segmented MRI slice"
                className="w-full h-auto"
              />
            )}
          </ModalBody>
          <ModalFooter>
            <Button onPress={() => setImageModalOpen(false)}>Close</Button>
          </ModalFooter>
        </ModalContent>
      </Modal>
    </div>
  );
}