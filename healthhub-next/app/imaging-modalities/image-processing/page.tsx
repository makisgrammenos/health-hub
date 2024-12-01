'use client';

import React, { useState, useRef, useEffect } from 'react';
import { Button } from '@nextui-org/button';
import { Card, CardBody } from '@nextui-org/card';
import { Slider } from '@nextui-org/slider';
import { Input } from '@nextui-org/input';
import { Image } from '@nextui-org/image';
import { Progress } from '@nextui-org/progress';
import debounce from 'lodash.debounce';

export default function MedicalImageEnhancement() {
  const [image, setImage] = useState<File | null>(null);
  const [uploadedImagePreview, setUploadedImagePreview] = useState<string | null>(null);
  const [processedSteps, setProcessedSteps] = useState<string[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const websocketRef = useRef<WebSocket | null>(null);

  // Parameters for preprocessing
  const [params, setParams] = useState({
    denoiseMethod: 'bilateral', // Denoising
    bilateralD: 9,
    bilateralSigmaColor: 75,
    bilateralSigmaSpace: 75,
    gaussianKernelSize: 5, // Gaussian Blur
    clipLimit: 2.0, // CLAHE
    tileGridSize: 8,
    unsharpStrength: 1.5, // Unsharp Masking
    normalizeMin: 0,
    normalizeMax: 255,
  });

  // WebSocket connection setup
  useEffect(() => {
    websocketRef.current = new WebSocket('ws://localhost:8000/imaging/image-processing/ws');
    websocketRef.current.onmessage = (event) => {
      if (event.data instanceof Blob) {
        const reader = new FileReader();
        reader.onload = () => {
          setProcessedSteps((prev) => [...prev, reader.result as string]);
          setLoading(false);
        };
        reader.readAsDataURL(event.data);
      } else {
        // Handle JSON messages
        try {
          const data = JSON.parse(event.data);
          if (data.message === 'Processing complete') {
            setLoading(false);
          }
        } catch (e) {
          console.error('Failed to parse message', e);
        }
      }
    };

    websocketRef.current.onclose = () => {
      console.log('WebSocket closed.');
    };

    websocketRef.current.onerror = (err) => {
      setError('An error occurred during WebSocket communication.');
      console.error(err);
    };

    return () => {
      websocketRef.current?.close();
    };
  }, []);

  // Handle image upload
  const handleImageUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const selectedFile = event.target.files?.[0];
    if (selectedFile) {
      setError(null);
      setImage(selectedFile);
      setUploadedImagePreview(URL.createObjectURL(selectedFile));
      setProcessedSteps([]);
    }
  };

  // Send image and parameters via WebSocket
  const startEnhancement = () => {
    if (!image) {
      setError('Please upload an image.');
      return;
    }

    setLoading(true);
    setError(null);
    setProcessedSteps([]);

    const reader = new FileReader();
    reader.onload = () => {
      const base64Image = reader.result?.toString().split(',')[1];
      websocketRef.current?.send(
        JSON.stringify({
          image: base64Image,
          params,
        })
      );
    };
    reader.readAsDataURL(image);
  };

  // Debounce the startEnhancement function to prevent excessive calls
  const debouncedEnhancement = useRef(debounce(startEnhancement, 500)).current;

  // Update parameter values and trigger enhancement
  const handleParamChange = (name: string, value: number | string) => {
    setParams((prev) => ({
      ...prev,
      [name]: value,
    }));
  };

  // Trigger enhancement when parameters change
  useEffect(() => {
    if (image) {
      debouncedEnhancement();
    }
    // Cleanup function to cancel debounce on unmount
    return () => {
      debouncedEnhancement.cancel();
    };
  }, [params]);

  return (
    <div className="max-w-7xl mx-auto p-8">
      <h1 className="text-4xl font-extrabold text-center mb-8">🩺 Medical Image Enhancement</h1>
      <p className="text-lg text-gray-600 text-center mb-12">
        Upload your medical image and adjust preprocessing parameters for real-time enhancement.
      </p>

      <section className="mb-12">
        <h2 className="text-2xl font-bold mb-6">Upload an Image</h2>
        <input
          type="file"
          accept="image/*"
          onChange={handleImageUpload}
          className="block w-full max-w-md mx-auto text-sm text-gray-500 file:mr-4 file:py-2 file:px-4 file:rounded-full file:border-0 file:bg-indigo-50 file:text-indigo-700 hover:file:bg-indigo-100 mb-4"
        />
        {uploadedImagePreview && (
          <Card className="mb-6">
            <CardBody>
              <Image src={uploadedImagePreview} alt="Uploaded Image" width="600" height="400" className="rounded-lg" />
            </CardBody>
          </Card>
        )}
      </section>

      {image && (
        <section className="mb-12">
          <h2 className="text-2xl font-bold mb-6">Adjust Parameters for Preprocessing Steps</h2>
          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-6">
            <div>
              <h3 className="font-semibold">Denoising</h3>
              <Input
                label="Method"
                value={params.denoiseMethod}
                onChange={(e) => handleParamChange('denoiseMethod', e.target.value)}
              />
              <Slider
                label="Bilateral Diameter"
                value={params.bilateralD}
                min={1}
                max={50}
                onChange={(value) => handleParamChange('bilateralD', value)}
              />
              <Slider
                label="Sigma Color"
                value={params.bilateralSigmaColor}
                min={10}
                max={150}
                onChange={(value) => handleParamChange('bilateralSigmaColor', value)}
              />
              <Slider
                label="Sigma Space"
                value={params.bilateralSigmaSpace}
                min={10}
                max={150}
                onChange={(value) => handleParamChange('bilateralSigmaSpace', value)}
              />
            </div>
            <div>
              <h3 className="font-semibold">Gaussian Blur</h3>
              <Slider
                label="Kernel Size"
                value={params.gaussianKernelSize}
                min={1}
                max={15}
                step={2}
                onChange={(value) => handleParamChange('gaussianKernelSize', value)}
              />
            </div>
            <div>
              <h3 className="font-semibold">CLAHE (Contrast Enhancement)</h3>
              <Slider
                label="Clip Limit"
                value={params.clipLimit}
                min={0.1}
                max={5.0}
                step={0.1}
                onChange={(value) => handleParamChange('clipLimit', value)}
              />
              <Slider
                label="Tile Grid Size"
                value={params.tileGridSize}
                min={1}
                max={16}
                onChange={(value) => handleParamChange('tileGridSize', value)}
              />
            </div>
            <div>
              <h3 className="font-semibold">Unsharp Masking</h3>
              <Slider
                label="Strength"
                value={params.unsharpStrength}
                min={0.1}
                max={3.0}
                step={0.1}
                onChange={(value) => handleParamChange('unsharpStrength', value)}
              />
            </div>
          </div>
        </section>
      )}

      {loading && (
        <div className="mt-6">
          <Progress indeterminate />
          <p className="text-center text-gray-600 mt-2">Processing image...</p>
        </div>
      )}

      {processedSteps.length > 0 && (
        <section className="mt-12">
          <h2 className="text-2xl font-bold mb-6">Enhanced Steps (Real-Time)</h2>
          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-6">
            {processedSteps.map((step, idx) => (
              <Card key={idx}>
                <CardBody>
                  <Image src={step} alt={`Step ${idx + 1}`} width="600" height="400" className="rounded-lg" />
                  <p className="mt-2 text-center text-sm">Step {idx + 1}</p>
                </CardBody>
              </Card>
            ))}
          </div>
        </section>
      )}

      {error && <p className="text-red-600 mt-4 text-center">{error}</p>}
    </div>
  );
}
