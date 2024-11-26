'use client';

import React, { useState } from 'react';
import axios from 'axios';
import { Card, CardBody, CardFooter } from '@nextui-org/card';
import { Button } from '@nextui-org/button';
import { Progress } from '@nextui-org/progress';

export default function BrainTumorSegmentation() {
  const [modalities, setModalities] = useState({
    T1c: null,
    T1: null,
    T2: null,
    FLAIR: null,
  });
  const [segmentationResult, setSegmentationResult] = useState<{
    segmented_slices: string[];
    download_url: string;
  } | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleFileChange = (modality: string, file: File) => {
    setModalities((prev) => ({ ...prev, [modality]: file }));
  };

  const startSegmentation = async () => {
    if (!modalities.T1c || !modalities.T1 || !modalities.T2 || !modalities.FLAIR) {
      setError('Please upload all four modalities.');
      return;
    }

    setLoading(true);
    setError(null);

    const formData = new FormData();
    Object.entries(modalities).forEach(([key, file]) => {
      formData.append(key, file as Blob);
    });

    try {
      const response = await axios.post('http://localhost:8000/imaging/brain-tumor/segment', formData, {
        headers: { 'Content-Type': 'multipart/form-data' },
      });
      setSegmentationResult(response.data);
    } catch (err) {
      setError('An error occurred during segmentation. Please try again.');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="max-w-7xl mx-auto p-8">
      <h1 className="text-4xl font-extrabold text-center mb-8">🧠 Brain Tumor Segmentation</h1>
      <p className="text-lg text-gray-600 text-center mb-12">
        Upload MRI modalities, perform segmentation, and view results step by step.
      </p>

      {/* Step 1: Upload Files */}
      <section className="mb-12">
        <h2 className="text-2xl font-bold mb-6">Step 1: Upload MRI Modalities</h2>
        <p className="mb-6 text-gray-600">
          Upload the four MRI modalities (T1c, T1, T2, FLAIR) as NIfTI (.nii, .nii.gz) files.
        </p>
        <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-4 gap-6">
          {['T1c', 'T1', 'T2', 'FLAIR'].map((modality) => (
            <Card key={modality} shadow="sm" className="flex flex-col items-center p-4">
              <CardBody>
                <h3 className="font-semibold text-center">{modality}</h3>
                <input
                  type="file"
                  accept=".nii,.nii.gz"
                  onChange={(e) => handleFileChange(modality, e.target.files?.[0] || null)}
                  className="block w-full text-sm text-gray-500 file:mr-4 file:py-2 file:px-4 file:rounded-full file:border-0 file:bg-indigo-50 file:text-indigo-700 hover:file:bg-indigo-100 mt-4"
                />
              </CardBody>
            </Card>
          ))}
        </div>
        <div className="text-center mt-8">
          <Button
            isDisabled={loading}
            isLoading={loading}
            onPress={startSegmentation}
            color="primary"
            size="lg"
          >
            Start Segmentation
          </Button>
        </div>
        {error && <div className="text-red-600 mt-4 text-center">{error}</div>}
      </section>

      {/* Step 2: Segmentation Progress */}
      <section className="mb-12">
        <h2 className="text-2xl font-bold mb-6">Step 2: Segmentation Progress</h2>
        {loading ? (
          <Progress color="primary" value={100} label="Processing... Please wait." />
        ) : (
          <p className="text-gray-600">No processing in progress.</p>
        )}
      </section>

      {/* Step 3: Segmentation Results */}
      <section>
        <h2 className="text-2xl font-bold mb-6">Step 3: Segmentation Results</h2>
        {segmentationResult ? (
          <div>
            <p className="text-gray-600 text-center mb-6">
              View the segmented slices below or download the full segmentation results.
            </p>
            <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-4 gap-6">
              {segmentationResult.segmented_slices.map((slice, idx) => (
                <Card key={idx} shadow="sm" className="flex flex-col items-center p-4">
                  <CardBody>
                    <img
                      src={slice}
                      alt={`Segmented Slice ${idx + 1}`}
                      className="w-full h-auto rounded-lg"
                    />
                  </CardBody>
                  <CardFooter>
                    <p className="text-center text-sm font-semibold">Slice {idx + 1}</p>
                  </CardFooter>
                </Card>
              ))}
            </div>
            <div className="text-center mt-8">
              <Button
                as="a"
                href={segmentationResult.download_url}
                download
                className="mt-4"
                color="success"
                size="lg"
              >
                Download Segmentation Results
              </Button>
            </div>
          </div>
        ) : (
          <p className="text-gray-600">No results available yet. Please complete segmentation.</p>
        )}
      </section>
    </div>
  );
}
