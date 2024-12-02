'use client';

import React, { useState, useRef, useEffect, useCallback, useMemo } from 'react';
import { Slider } from '@nextui-org/slider';
import { Progress } from '@nextui-org/progress';
import { Switch } from '@nextui-org/switch';
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
  });

  // Methods to apply
  const [methods, setMethods] = useState({
    denoising: true,
    normalization: true,
    histogramEqualization: true,
    clahe: true,
    unsharpMasking: true,
  });

  // WebSocket connection setup
  useEffect(() => {
    websocketRef.current = new WebSocket('ws://localhost:8000/imaging/image-processing/ws');

    websocketRef.current.onopen = () => console.log('WebSocket connected.');
    websocketRef.current.onmessage = (event) => {
      if (event.data instanceof Blob) {
        // Handle Blob data from the WebSocket
        const reader = new FileReader();
        reader.onload = () => {
          const base64Image = reader.result as string;
          console.log('Received Base64 Image:', base64Image); // Debugging: Log received image
          setProcessedImage(base64Image);
          setLoading(false);
        };
        reader.readAsDataURL(event.data);
      } else if (typeof event.data === 'string') {
        // Handle JSON messages
        try {
          const data = JSON.parse(event.data);
          if (data.message === 'Processing complete') {
            setLoading(false);
          } else if (data.error) {
            setError(data.error);
            setLoading(false);
          }
        } catch (e) {
          console.error('Failed to parse message:', e);
        }
      }
    };

    websocketRef.current.onclose = () => console.log('WebSocket closed.');
    websocketRef.current.onerror = (err) => {
      setError('An error occurred during WebSocket communication.');
      console.error('WebSocket error:', err);
    };

    return () => websocketRef.current?.close();
  }, []);

  const handleImageUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const selectedFile = event.target.files?.[0];
    if (selectedFile) {
      setError(null);
      setImage(selectedFile);
      setUploadedImagePreview(URL.createObjectURL(selectedFile));
      setProcessedImage(null);
    }
  };

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
    setError(null);
    setProcessedImage(null);

    const reader = new FileReader();
    reader.onload = () => {
      const base64Image = reader.result?.toString().split(',')[1];
      websocketRef.current?.send(
        JSON.stringify({
          image: base64Image,
          params,
          methods,
        })
      );
    };
    reader.readAsDataURL(image);
  }, [image, params, methods]);

  const debouncedEnhancement = useMemo(() => debounce(startEnhancement, 500), [startEnhancement]);

  const handleParamChange = (name: string, value: number | string) => {
    setParams((prev) => ({
      ...prev,
      [name]: value,
    }));
  };

  const handleMethodChange = (name: string, value: boolean) => {
    setMethods((prev) => ({
      ...prev,
      [name]: value,
    }));
  };

  useEffect(() => {
    if (image) {
      debouncedEnhancement();
    }
    return () => {
      debouncedEnhancement.cancel();
    };
  }, [params, methods, debouncedEnhancement, image]);

  return (
    <div className="flex h-screen bg-gray-100">
      {/* Sidebar */}
      <aside className="w-1/4 bg-gray-800 text-white p-6 overflow-y-auto shadow-md">
        <h2 className="text-lg font-bold mb-4">Preprocessing Methods</h2>
        {Object.entries(methods).map(([method, isActive]) => (
          <div key={method} className="flex items-center justify-between mb-4">
            <span className="capitalize">{method.replace(/([A-Z])/g, ' $1')}</span>
            <Switch isSelected={isActive} onChange={(e) => handleMethodChange(method, e.target.checked)} />
          </div>
        ))}
      </aside>

      {/* Main Panel */}
      <main className="flex-1 flex flex-col">
        {/* Top Section */}
        <div className="p-6 bg-white shadow-md border-b">
          <h1 className="text-2xl font-bold mb-4">🩺 Medical Image Enhancement</h1>
          <div className="flex items-center gap-4 mb-6">
            <input
              type="file"
              accept="image/*"
              onChange={handleImageUpload}
              className="block w-full max-w-md text-sm text-gray-500 file:mr-4 file:py-2 file:px-4 file:rounded-full file:border-0 file:bg-blue-50 file:text-blue-700 hover:file:bg-blue-100"
            />
          </div>
          {uploadedImagePreview && (
            <div className="mb-6">
              <h2 className="font-semibold text-lg">Uploaded Image:</h2>
              <img
                src={uploadedImagePreview}
                alt="Uploaded Preview"
                className="w-full max-w-lg mx-auto rounded-md shadow-md"
              />
            </div>
          )}
        </div>

        {/* Bottom Section */}
        <div className="flex-1 p-6 overflow-y-auto bg-gray-50">
          {loading && (
            <div className="flex justify-center items-center h-full">
              <Progress indeterminate />
              <p className="text-gray-600 mt-2">Processing image...</p>
            </div>
          )}
          {processedImage && (
            <div>
              <h2 className="font-semibold text-lg mb-4">Enhanced Image:</h2>
              <img
                src={processedImage}
                alt="Processed Image"
                className="w-full max-w-lg mx-auto rounded-md shadow-md"
              />
            </div>
          )}
          {error && <p className="text-red-600 mt-4">{error}</p>}
        </div>
      </main>
    </div>
  );
}
