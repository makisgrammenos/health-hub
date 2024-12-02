'use client';

import React, { useState, useRef, useEffect, useCallback, useMemo } from 'react';
import { Card, CardHeader, CardBody } from '@nextui-org/card';
import { Input } from '@nextui-org/input';
import { Progress } from '@nextui-org/progress';
import { Slider } from '@nextui-org/slider';
import { Checkbox, CheckboxGroup } from '@nextui-org/checkbox'; // Checkbox for toggling methods
import { Select, SelectItem } from '@nextui-org/select'; // Dropdown component
import debounce from 'lodash.debounce';

export default function MedicalImageEnhancement() {
  const [image, setImage] = useState<File | null>(null);
  const [uploadedImagePreview, setUploadedImagePreview] = useState<string | null>(null);
  const [processedImage, setProcessedImage] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const websocketRef = useRef<WebSocket | null>(null);

  // Parameters for preprocessing
  const [params, setParams] = useState({
    denoiseMethod: 'bilateral',
    bilateralD: 9,
    bilateralSigmaColor: 75,
    bilateralSigmaSpace: 75,
    gaussianKernelSize: 5,
    clipLimit: 2.0,
    tileGridSize: 8,
    unsharpStrength: 1.5,
    normalizeMin: 0,
    normalizeMax: 255,
    windowWidth: 255,
    windowCenter: 127,
    thresholdMin: 50,
    thresholdMax: 200,
    pseudocolorMap: 'COLORMAP_JET',
    kernelSize: 3,
    morphOperation: 'dilation',
  });

  // Active preprocessing methods
  const [activeMethods, setActiveMethods] = useState<string[]>([]);

  // WebSocket connection setup
  useEffect(() => {
    websocketRef.current = new WebSocket('ws://localhost:8000/ws');

    websocketRef.current.onopen = () => console.log('WebSocket connected.');

    websocketRef.current.onmessage = (event) => {
      if (typeof event.data === 'string') {
        setProcessedImage(`data:image/png;base64,${event.data}`);
        setLoading(false);
      }
    };

    websocketRef.current.onerror = (err) => {
      setError('WebSocket error: ' + err);
      console.error(err);
    };

    websocketRef.current.onclose = () => console.log('WebSocket closed.');

    return () => websocketRef.current?.close();
  }, []);

  // Handle image upload
  const handleImageUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const selectedFile = event.target.files?.[0];
    if (selectedFile) {
      setImage(selectedFile);
      setUploadedImagePreview(URL.createObjectURL(selectedFile));
      setProcessedImage(null);
      setError(null);
    }
  };

  // Start enhancement process
  const startEnhancement = useCallback(() => {
    if (!image) {
      setError('Please upload an image.');
      return;
    }
    if (!websocketRef.current || websocketRef.current.readyState !== WebSocket.OPEN) {
      setError('WebSocket connection is not open.');
      return;
    }
    setLoading(true);
    const reader = new FileReader();
    reader.onload = () => {
      websocketRef.current?.send(
        JSON.stringify({
          image: reader.result?.toString().split(',')[1],
          params,
          methods: activeMethods.reduce((acc, method) => {
            acc[method] = true;
            return acc;
          }, {}),
        })
      );
    };
    reader.readAsDataURL(image);
  }, [image, params, activeMethods]);

  // Debounce the enhancement function
  const debouncedEnhancement = useMemo(() => debounce(startEnhancement, 500), [startEnhancement]);

  // Handle parameter changes
  const handleParamChange = (name: string, value: any) => {
    setParams((prev) => ({
      ...prev,
      [name]: value,
    }));
    debouncedEnhancement(); // Trigger the enhancement process when parameters change
  };

  // Handle active methods change
  const handleMethodChange = (methods: string[]) => {
    setActiveMethods(methods);
    debouncedEnhancement(); // Trigger the enhancement process when methods change
  };

  return (
    <div className="flex flex-col h-screen">
      {/* Top Bar */}
      <div className="bg-gray-100 p-4 border-b">
        <h1 className="text-2xl font-bold text-center">🩺 Medical Image Enhancement</h1>
      </div>

      {/* Main Content */}
      <div className="flex flex-grow">
        {/* Sidebar for Controls */}
        <aside className="w-1/4 bg-gray-50 p-4 border-r overflow-y-auto">
          <Input
            type="file"
            accept="image/*"
            onChange={handleImageUpload}
            label="Upload Image"
            fullWidth
            className="mb-6"
          />

          <h3 className="text-lg font-semibold mb-4">Preprocessing Methods</h3>
          <CheckboxGroup
            value={activeMethods}
            onChange={(values) => handleMethodChange(values)}
          >
            <Checkbox value="denoising">Denoising</Checkbox>
            <Checkbox value="normalization">Normalization</Checkbox>
            <Checkbox value="histogramEqualization">Histogram Equalization</Checkbox>
            <Checkbox value="clahe">CLAHE</Checkbox>
            <Checkbox value="unsharpMasking">Unsharp Masking</Checkbox>
            <Checkbox value="windowing">Windowing</Checkbox>
            <Checkbox value="thresholdSegmentation">Threshold Segmentation</Checkbox>
            <Checkbox value="pseudocolor">Pseudocolor</Checkbox>
            <Checkbox value="morphologicalOperations">Morphological Operations</Checkbox>
          </CheckboxGroup>

          {/* Adjustable Parameters */}
          <div className="mt-6">
            {activeMethods.includes('denoising') && (
              <>
                <Select
                  label="Denoise Method"
                  value={params.denoiseMethod}
                  onChange={(value) => handleParamChange('denoiseMethod', value)}
                >
                  <SelectItem key="bilateral" value="bilateral">Bilateral</SelectItem>
                  <SelectItem key="gaussian" value="gaussian">Gaussian</SelectItem>
                  <SelectItem key="median" value="median">Median</SelectItem>
                </Select>
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
              </>
            )}
            {activeMethods.includes('clahe') && (
              <>
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
              </>
            )}
            {/* Add additional adjustable parameters for other methods here */}
          </div>
        </aside>

        {/* Main Image Area */}
        <main className="flex-1 p-4 flex flex-col items-center justify-center bg-gray-100">
          {loading && (
            <div className="flex flex-col items-center">
              <Progress indeterminate />
              <p className="text-gray-600 mt-2">Processing image...</p>
            </div>
          )}
          {!loading && processedImage && (
            <img
              src={processedImage}
              alt="Processed Image"
              className="w-full max-w-2xl rounded-lg shadow-md"
            />
          )}
          {!loading && !processedImage && uploadedImagePreview && (
            <img
              src={uploadedImagePreview}
              alt="Uploaded Image"
              className="w-full max-w-2xl rounded-lg shadow-md"
            />
          )}
          {error && <p className="text-red-500 mt-4">{error}</p>}
        </main>
      </div>
    </div>
  );
}
