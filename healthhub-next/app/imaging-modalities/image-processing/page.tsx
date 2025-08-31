'use client';

import React, { useState, useRef, useEffect, useCallback } from 'react';

/**
 * Medical Image Enhancement Component
 *
 * A comprehensive interface for applying various image processing techniques
 * to enhance medical imaging data. This component enables researchers and
 * clinicians to apply real-time adjustments and compare preprocessing methods
 * for improved diagnostic visualization.
 *
 * @returns {JSX.Element} The Medical Image Enhancement interface
 */
export default function MedicalImageEnhancement(): JSX.Element {
  // Image state management
  const [image, setImage] = useState<File | null>(null);
  const [uploadedImagePreview, setUploadedImagePreview] = useState<string | null>(null);
  const [processedImage, setProcessedImage] = useState<string | null>(null);
  const [originalImage, setOriginalImage] = useState<string | null>(null);
  const [loading, setLoading] = useState<boolean>(false);
  const [error, setError] = useState<string | null>(null);
  const websocketRef = useRef<WebSocket | null>(null);

  // Parameters for image processing techniques with appropriate defaults
  const [params, setParams] = useState({
    // Denoising parameters
    denoiseMethod: 'bilateral',
    bilateralD: 9,
    bilateralSigmaColor: 75,
    bilateralSigmaSpace: 75,
    gaussianKernelSize: 5,
    medianKernelSize: 5,

    // CLAHE parameters
    clipLimit: 2.0,
    tileGridSize: 8,

    // Sharpening parameters
    unsharpStrength: 1.5,

    // Intensity transformation parameters
    normalizeMin: 0,
    normalizeMax: 255,

    // Windowing parameters (for radiological images)
    windowWidth: 255,
    windowCenter: 127,

    // Segmentation parameters
    thresholdMin: 50,
    thresholdMax: 200,

    // Visualization parameters
    pseudocolorMap: 'COLORMAP_JET',

    // Morphological parameters
    kernelSize: 3,
    morphOperation: 'dilation',
  });

  // Available preprocessing methods with descriptions
  const methodDescriptions: Record<string, string> = {
    denoising: "Reduces noise in medical images to improve signal-to-noise ratio",
    normalization: "Standardizes pixel intensity range for consistent visualization",
    histogramEqualization: "Enhances contrast by equalizing the image histogram",
    clahe: "Contrast Limited Adaptive Histogram Equalization for local contrast enhancement",
    unsharpMasking: "Enhances edges and fine details by subtracting a blurred version",
    windowing: "Adjusts visibility of structures in different radiodensity ranges",
    thresholdSegmentation: "Segments regions of interest based on intensity thresholds",
    pseudocolor: "Applies color mapping to grayscale images for enhanced visualization",
    morphologicalOperations: "Applies structural transformations for feature extraction or noise removal",
  };

  // Track active preprocessing methods
  const [activeMethods, setActiveMethods] = useState<string[]>([]);

  // WebSocket connection setup
  useEffect(() => {
    try {
      websocketRef.current = new WebSocket('ws://localhost:8000/imaging/image-processing/ws');

      websocketRef.current.onopen = () => {
        console.log('WebSocket connection established successfully.');
      };

      websocketRef.current.onmessage = (event) => {
        if (typeof event.data === 'string') {
          try {
            const data = JSON.parse(event.data);
            if (data.error) {
              setError(`Processing error: ${data.error}`);
              setLoading(false);
            }
          } catch (e) {
            setProcessedImage(`data:image/png;base64,${event.data}`);
            setLoading(false);
          }
        }
      };

      websocketRef.current.onerror = (err) => {
        console.error('WebSocket error:', err);
        setLoading(false);
      };

      websocketRef.current.onclose = () => {
        console.log('WebSocket connection closed.');
      };
    } catch (err) {
      setError('Failed to establish WebSocket connection.');
      console.error('WebSocket initialization error:', err);
    }

    return () => {
      if (websocketRef.current) {
        websocketRef.current.close();
      }
    };
  }, []);

  const startEnhancement = useCallback(() => {
    if (!image) {
      setError('Please upload a medical image to process.');
      return;
    }

    if (!websocketRef.current || websocketRef.current.readyState !== WebSocket.OPEN) {
      setError('WebSocket connection is not established. Please refresh the page.');
      return;
    }

    if (activeMethods.length === 0) {
      setProcessedImage(originalImage);
      return;
    }

    setLoading(true);
    setError(null);

    const reader = new FileReader();
    reader.onload = () => {
      if (reader.result && typeof reader.result === 'string' && websocketRef.current) {
        const base64Data = reader.result.split(',')[1];

        websocketRef.current.send(
          JSON.stringify({
            image: base64Data,
            params,
            methods: activeMethods.reduce((acc, method) => {
              acc[method] = true;
              return acc;
            }, {} as Record<string, boolean>),
          })
        );
      } else {
        setError('Error reading image data.');
        setLoading(false);
      }
    };

    reader.onerror = () => {
      setError('Failed to read the image file.');
      setLoading(false);
    };

    reader.readAsDataURL(image);
  }, [image, params, activeMethods, originalImage]);

  useEffect(() => {
    if (image) {
      const debounceTimeout = setTimeout(() => {
        startEnhancement();
      }, 500);

      return () => clearTimeout(debounceTimeout);
    }
  }, [image, params, activeMethods, startEnhancement]);

  const handleImageUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const selectedFile = event.target.files?.[0];

    if (selectedFile) {
      const validImageTypes = ['image/jpeg', 'image/png', 'image/dicom', 'image/tiff'];
      if (!validImageTypes.includes(selectedFile.type) && !selectedFile.name.endsWith('.dcm')) {
        setError('Please upload a valid medical image format (JPEG, PNG, DICOM, TIFF).');
        return;
      }

      const objectUrl = URL.createObjectURL(selectedFile);
      setImage(selectedFile);
      setUploadedImagePreview(objectUrl);
      setOriginalImage(objectUrl);
      setProcessedImage(null);
      setError(null);
      setActiveMethods([]);
    }
  };

  const handleParamChange = (name: string, value: any) => {
    setParams((prev) => ({
      ...prev,
      [name]: value,
    }));
  };

  const handleMethodChange = (method: string, isChecked: boolean) => {
    setActiveMethods((prev) => {
      return isChecked
        ? [...prev, method]
        : prev.filter((m) => m !== method);
    });
  };

  const handleResetParams = () => {
    setParams({
      denoiseMethod: 'bilateral',
      bilateralD: 9,
      bilateralSigmaColor: 75,
      bilateralSigmaSpace: 75,
      gaussianKernelSize: 5,
      medianKernelSize: 5,
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
  };

  const toggleAllMethods = (enable: boolean) => {
    if (enable) {
      setActiveMethods(Object.keys(methodDescriptions));
    } else {
      setActiveMethods([]);
    }
  };

  return (
    <div className="flex flex-col h-screen bg-gray-50">
      <header className="p-4 bg-blue-800 text-white">
        <h1 className="text-2xl font-bold">Medical Image Enhancement and Analysis</h1>
        <p className="text-sm">Advanced visualization and preprocessing tools for clinical and research applications</p>
      </header>

      <div className="flex flex-grow overflow-hidden">
        {/* Sidebar for Controls */}
        <aside className="w-1/3 p-4 border-r overflow-y-auto bg-white">
          <div className="bg-white rounded-lg shadow-md mb-4 p-4">
            <h2 className="text-lg font-bold mb-4">Image Input</h2>
            <div className="w-full">
              <label className="block text-sm font-medium text-gray-700 mb-2">
                Upload Medical Image
              </label>
              <input
                type="file"
                accept="image/*,.dcm"
                onChange={handleImageUpload}
                className="block w-full text-sm text-gray-900 border border-gray-300 rounded-lg cursor-pointer bg-gray-50 focus:outline-none focus:border-blue-500 file:mr-4 file:py-2 file:px-4 file:rounded-md file:border-0 file:text-sm file:font-semibold file:bg-blue-50 file:text-blue-700 hover:file:bg-blue-100"
              />
              <p className="mt-1 text-sm text-gray-500">
                Supported formats: JPEG, PNG, DICOM, TIFF
              </p>
            </div>

            <div className="flex space-x-2 mt-4">
              <button
                onClick={() => toggleAllMethods(true)}
                disabled={!image}
                className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:bg-gray-300 disabled:cursor-not-allowed text-sm"
              >
                Enable All
              </button>
              <button
                onClick={() => toggleAllMethods(false)}
                disabled={!image}
                className="px-4 py-2 bg-red-600 text-white rounded-md hover:bg-red-700 disabled:bg-gray-300 disabled:cursor-not-allowed text-sm"
              >
                Disable All
              </button>
              <button
                onClick={handleResetParams}
                disabled={!image}
                className="px-4 py-2 bg-yellow-600 text-white rounded-md hover:bg-yellow-700 disabled:bg-gray-300 disabled:cursor-not-allowed text-sm"
              >
                Reset Params
              </button>
            </div>
          </div>

          <div className="bg-white rounded-lg shadow-md p-4">
            <h2 className="text-lg font-bold mb-2">Processing Methods</h2>
            <p className="text-sm text-gray-500 mb-4">Select techniques to enhance the image</p>

            {/* Denoising */}
            <div className="mb-4 pb-4 border-b">
              <div className="flex items-center mb-2">
                <input
                  type="checkbox"
                  id="denoising"
                  checked={activeMethods.includes('denoising')}
                  onChange={(e) => handleMethodChange('denoising', e.target.checked)}
                  className="mr-2"
                />
                <label htmlFor="denoising" className="font-medium cursor-pointer">Denoising</label>
                <span className="ml-2 text-gray-500 cursor-help text-sm" title={methodDescriptions.denoising}>ⓘ</span>
              </div>

              {activeMethods.includes('denoising') && (
                <div className="pl-6 mt-2 space-y-3">
                  <div>
                    <label className="block text-sm font-medium text-gray-700 mb-1">Filter Algorithm</label>
                    <select
                      value={params.denoiseMethod}
                      onChange={(e) => handleParamChange('denoiseMethod', e.target.value)}
                      className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                    >
                      <option value="bilateral">Bilateral Filter</option>
                      <option value="gaussian">Gaussian Filter</option>
                      <option value="median">Median Filter</option>
                    </select>
                  </div>

                  {params.denoiseMethod === 'bilateral' && (
                    <>
                      <div>
                        <label className="block text-sm text-gray-700 mb-1">Diameter: {params.bilateralD}</label>
                        <input
                          type="range"
                          min="3"
                          max="25"
                          step="2"
                          value={params.bilateralD}
                          onChange={(e) => handleParamChange('bilateralD', Number(e.target.value))}
                          className="w-full"
                        />
                      </div>
                      <div>
                        <label className="block text-sm text-gray-700 mb-1">Color Sigma: {params.bilateralSigmaColor}</label>
                        <input
                          type="range"
                          min="10"
                          max="150"
                          step="5"
                          value={params.bilateralSigmaColor}
                          onChange={(e) => handleParamChange('bilateralSigmaColor', Number(e.target.value))}
                          className="w-full"
                        />
                      </div>
                      <div>
                        <label className="block text-sm text-gray-700 mb-1">Space Sigma: {params.bilateralSigmaSpace}</label>
                        <input
                          type="range"
                          min="10"
                          max="150"
                          step="5"
                          value={params.bilateralSigmaSpace}
                          onChange={(e) => handleParamChange('bilateralSigmaSpace', Number(e.target.value))}
                          className="w-full"
                        />
                      </div>
                    </>
                  )}

                  {params.denoiseMethod === 'gaussian' && (
                    <div>
                      <label className="block text-sm text-gray-700 mb-1">Kernel Size: {params.gaussianKernelSize}</label>
                      <input
                        type="range"
                        min="1"
                        max="15"
                        step="2"
                        value={params.gaussianKernelSize}
                        onChange={(e) => handleParamChange('gaussianKernelSize', Number(e.target.value))}
                        className="w-full"
                      />
                    </div>
                  )}

                  {params.denoiseMethod === 'median' && (
                    <div>
                      <label className="block text-sm text-gray-700 mb-1">Kernel Size: {params.medianKernelSize}</label>
                      <input
                        type="range"
                        min="1"
                        max="15"
                        step="2"
                        value={params.medianKernelSize}
                        onChange={(e) => handleParamChange('medianKernelSize', Number(e.target.value))}
                        className="w-full"
                      />
                    </div>
                  )}
                </div>
              )}
            </div>

            {/* Normalization */}
            <div className="mb-4 pb-4 border-b">
              <div className="flex items-center mb-2">
                <input
                  type="checkbox"
                  id="normalization"
                  checked={activeMethods.includes('normalization')}
                  onChange={(e) => handleMethodChange('normalization', e.target.checked)}
                  className="mr-2"
                />
                <label htmlFor="normalization" className="font-medium cursor-pointer">Intensity Normalization</label>
                <span className="ml-2 text-gray-500 cursor-help text-sm" title={methodDescriptions.normalization}>ⓘ</span>
              </div>

              {activeMethods.includes('normalization') && (
                <div className="pl-6 mt-2 space-y-3">
                  <div>
                    <label className="block text-sm text-gray-700 mb-1">Min Intensity: {params.normalizeMin}</label>
                    <input
                      type="range"
                      min="0"
                      max="255"
                      value={params.normalizeMin}
                      onChange={(e) => handleParamChange('normalizeMin', Number(e.target.value))}
                      className="w-full"
                    />
                  </div>
                  <div>
                    <label className="block text-sm text-gray-700 mb-1">Max Intensity: {params.normalizeMax}</label>
                    <input
                      type="range"
                      min="0"
                      max="255"
                      value={params.normalizeMax}
                      onChange={(e) => handleParamChange('normalizeMax', Number(e.target.value))}
                      className="w-full"
                    />
                  </div>
                </div>
              )}
            </div>

            {/* Histogram Equalization */}
            <div className="mb-4 pb-4 border-b">
              <div className="flex items-center mb-2">
                <input
                  type="checkbox"
                  id="histogramEqualization"
                  checked={activeMethods.includes('histogramEqualization')}
                  onChange={(e) => handleMethodChange('histogramEqualization', e.target.checked)}
                  className="mr-2"
                />
                <label htmlFor="histogramEqualization" className="font-medium cursor-pointer">Histogram Equalization</label>
                <span className="ml-2 text-gray-500 cursor-help text-sm" title={methodDescriptions.histogramEqualization}>ⓘ</span>
              </div>

              {activeMethods.includes('histogramEqualization') && (
                <div className="pl-6 mt-2">
                  <p className="text-sm text-gray-600">Global histogram equalization with no adjustable parameters.</p>
                </div>
              )}
            </div>

            {/* CLAHE */}
            <div className="mb-4 pb-4 border-b">
              <div className="flex items-center mb-2">
                <input
                  type="checkbox"
                  id="clahe"
                  checked={activeMethods.includes('clahe')}
                  onChange={(e) => handleMethodChange('clahe', e.target.checked)}
                  className="mr-2"
                />
                <label htmlFor="clahe" className="font-medium cursor-pointer">CLAHE</label>
                <span className="ml-2 text-gray-500 cursor-help text-sm" title={methodDescriptions.clahe}>ⓘ</span>
              </div>

              {activeMethods.includes('clahe') && (
                <div className="pl-6 mt-2 space-y-3">
                  <div>
                    <label className="block text-sm text-gray-700 mb-1">Clip Limit: {params.clipLimit.toFixed(1)}</label>
                    <input
                      type="range"
                      min="0.5"
                      max="8.0"
                      step="0.1"
                      value={params.clipLimit}
                      onChange={(e) => handleParamChange('clipLimit', Number(e.target.value))}
                      className="w-full"
                    />
                  </div>
                  <div>
                    <label className="block text-sm text-gray-700 mb-1">Tile Grid Size: {params.tileGridSize}</label>
                    <input
                      type="range"
                      min="2"
                      max="16"
                      step="2"
                      value={params.tileGridSize}
                      onChange={(e) => handleParamChange('tileGridSize', Number(e.target.value))}
                      className="w-full"
                    />
                  </div>
                </div>
              )}
            </div>

            {/* Edge Enhancement */}
            <div className="mb-4 pb-4 border-b">
              <div className="flex items-center mb-2">
                <input
                  type="checkbox"
                  id="unsharpMasking"
                  checked={activeMethods.includes('unsharpMasking')}
                  onChange={(e) => handleMethodChange('unsharpMasking', e.target.checked)}
                  className="mr-2"
                />
                <label htmlFor="unsharpMasking" className="font-medium cursor-pointer">Edge Enhancement</label>
                <span className="ml-2 text-gray-500 cursor-help text-sm" title={methodDescriptions.unsharpMasking}>ⓘ</span>
              </div>

              {activeMethods.includes('unsharpMasking') && (
                <div className="pl-6 mt-2 space-y-3">
                  <div>
                    <label className="block text-sm text-gray-700 mb-1">Strength: {params.unsharpStrength.toFixed(1)}</label>
                    <input
                      type="range"
                      min="0.1"
                      max="5.0"
                      step="0.1"
                      value={params.unsharpStrength}
                      onChange={(e) => handleParamChange('unsharpStrength', Number(e.target.value))}
                      className="w-full"
                    />
                  </div>
                </div>
              )}
            </div>

            {/* Windowing */}
            <div className="mb-4 pb-4 border-b">
              <div className="flex items-center mb-2">
                <input
                  type="checkbox"
                  id="windowing"
                  checked={activeMethods.includes('windowing')}
                  onChange={(e) => handleMethodChange('windowing', e.target.checked)}
                  className="mr-2"
                />
                <label htmlFor="windowing" className="font-medium cursor-pointer">Radiological Windowing</label>
                <span className="ml-2 text-gray-500 cursor-help text-sm" title={methodDescriptions.windowing}>ⓘ</span>
              </div>

              {activeMethods.includes('windowing') && (
                <div className="pl-6 mt-2 space-y-3">
                  <div>
                    <label className="block text-sm text-gray-700 mb-1">Window Width: {params.windowWidth}</label>
                    <input
                      type="range"
                      min="1"
                      max="512"
                      value={params.windowWidth}
                      onChange={(e) => handleParamChange('windowWidth', Number(e.target.value))}
                      className="w-full"
                    />
                  </div>
                  <div>
                    <label className="block text-sm text-gray-700 mb-1">Window Center: {params.windowCenter}</label>
                    <input
                      type="range"
                      min="1"
                      max="512"
                      value={params.windowCenter}
                      onChange={(e) => handleParamChange('windowCenter', Number(e.target.value))}
                      className="w-full"
                    />
                  </div>
                </div>
              )}
            </div>

            {/* Threshold Segmentation */}
            <div className="mb-4 pb-4 border-b">
              <div className="flex items-center mb-2">
                <input
                  type="checkbox"
                  id="thresholdSegmentation"
                  checked={activeMethods.includes('thresholdSegmentation')}
                  onChange={(e) => handleMethodChange('thresholdSegmentation', e.target.checked)}
                  className="mr-2"
                />
                <label htmlFor="thresholdSegmentation" className="font-medium cursor-pointer">Threshold Segmentation</label>
                <span className="ml-2 text-gray-500 cursor-help text-sm" title={methodDescriptions.thresholdSegmentation}>ⓘ</span>
              </div>

              {activeMethods.includes('thresholdSegmentation') && (
                <div className="pl-6 mt-2 space-y-3">
                  <div>
                    <label className="block text-sm text-gray-700 mb-1">Lower Threshold: {params.thresholdMin}</label>
                    <input
                      type="range"
                      min="0"
                      max="255"
                      value={params.thresholdMin}
                      onChange={(e) => handleParamChange('thresholdMin', Number(e.target.value))}
                      className="w-full"
                    />
                  </div>
                  <div>
                    <label className="block text-sm text-gray-700 mb-1">Upper Threshold: {params.thresholdMax}</label>
                    <input
                      type="range"
                      min="0"
                      max="255"
                      value={params.thresholdMax}
                      onChange={(e) => handleParamChange('thresholdMax', Number(e.target.value))}
                      className="w-full"
                    />
                  </div>
                </div>
              )}
            </div>

            {/* Pseudocolor */}
            <div className="mb-4 pb-4 border-b">
              <div className="flex items-center mb-2">
                <input
                  type="checkbox"
                  id="pseudocolor"
                  checked={activeMethods.includes('pseudocolor')}
                  onChange={(e) => handleMethodChange('pseudocolor', e.target.checked)}
                  className="mr-2"
                />
                <label htmlFor="pseudocolor" className="font-medium cursor-pointer">Pseudocolor Mapping</label>
                <span className="ml-2 text-gray-500 cursor-help text-sm" title={methodDescriptions.pseudocolor}>ⓘ</span>
              </div>

              {activeMethods.includes('pseudocolor') && (
                <div className="pl-6 mt-2">
                  <label className="block text-sm font-medium text-gray-700 mb-1">Color Map</label>
                  <select
                    value={params.pseudocolorMap}
                    onChange={(e) => handleParamChange('pseudocolorMap', e.target.value)}
                    className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                  >
                    <option value="COLORMAP_JET">Jet</option>
                    <option value="COLORMAP_VIRIDIS">Viridis</option>
                    <option value="COLORMAP_HOT">Hot</option>
                    <option value="COLORMAP_BONE">Bone</option>
                    <option value="COLORMAP_PLASMA">Plasma</option>
                    <option value="COLORMAP_INFERNO">Inferno</option>
                    <option value="COLORMAP_MAGMA">Magma</option>
                    <option value="COLORMAP_CIVIDIS">Cividis</option>
                    <option value="COLORMAP_RAINBOW">Rainbow</option>
                    <option value="COLORMAP_OCEAN">Ocean</option>
                    <option value="COLORMAP_TURBO">Turbo</option>
                  </select>
                </div>
              )}
            </div>

            {/* Morphological Operations */}
            <div className="mb-4">
              <div className="flex items-center mb-2">
                <input
                  type="checkbox"
                  id="morphologicalOperations"
                  checked={activeMethods.includes('morphologicalOperations')}
                  onChange={(e) => handleMethodChange('morphologicalOperations', e.target.checked)}
                  className="mr-2"
                />
                <label htmlFor="morphologicalOperations" className="font-medium cursor-pointer">Morphological Operations</label>
                <span className="ml-2 text-gray-500 cursor-help text-sm" title={methodDescriptions.morphologicalOperations}>ⓘ</span>
              </div>

              {activeMethods.includes('morphologicalOperations') && (
                <div className="pl-6 mt-2 space-y-3">
                  <div>
                    <label className="block text-sm font-medium text-gray-700 mb-1">Operation Type</label>
                    <select
                      value={params.morphOperation}
                      onChange={(e) => handleParamChange('morphOperation', e.target.value)}
                      className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                    >
                      <option value="dilation">Dilation</option>
                      <option value="erosion">Erosion</option>
                      <option value="opening">Opening</option>
                      <option value="closing">Closing</option>
                      <option value="gradient">Morphological Gradient</option>
                      <option value="tophat">Top Hat</option>
                      <option value="blackhat">Black Hat</option>
                    </select>
                  </div>
                  <div>
                    <label className="block text-sm text-gray-700 mb-1">Kernel Size: {params.kernelSize}</label>
                    <input
                      type="range"
                      min="1"
                      max="15"
                      step="2"
                      value={params.kernelSize}
                      onChange={(e) => handleParamChange('kernelSize', Number(e.target.value))}
                      className="w-full"
                    />
                  </div>
                </div>
              )}
            </div>
          </div>
        </aside>

        {/* Main Image Display Area */}
        <main className="flex-1 p-4 flex flex-col items-center justify-center bg-gray-100">
          {/* Processing indicator */}
          {loading && (
            <div className="w-full max-w-md p-4 bg-white rounded-lg shadow-md">
              <div className="flex flex-col items-center">
                <div className="w-full bg-gray-200 rounded-full h-2.5 mb-4">
                  <div className="bg-blue-600 h-2.5 rounded-full animate-pulse" style={{ width: '45%' }}></div>
                </div>
                <p className="text-center text-gray-700">
                  Applying image processing algorithms...
                </p>
              </div>
            </div>
          )}

          {/* Processed image display */}
          {!loading && processedImage && (
            <div className="w-full max-w-4xl bg-white rounded-lg shadow-md">
              <div className="p-4 border-b flex justify-between items-start">
                <div>
                  <h2 className="text-lg font-bold">Enhanced Image</h2>
                  <p className="text-sm text-gray-600">
                    Applied methods: {activeMethods.length > 0 ?
                      activeMethods.map(method => method.charAt(0).toUpperCase() + method.slice(1)).join(', ') :
                      'None (showing original)'}
                  </p>
                </div>
                <button
                  onClick={() => {
                    const link = document.createElement('a');
                    link.href = processedImage;
                    link.download = `enhanced_image_${Date.now()}.png`;
                    document.body.appendChild(link);
                    link.click();
                    document.body.removeChild(link);
                  }}
                  className="px-4 py-2 bg-green-600 text-white rounded-md hover:bg-green-700 flex items-center gap-2"
                >
                  <svg
                    xmlns="http://www.w3.org/2000/svg"
                    width="20"
                    height="20"
                    viewBox="0 0 24 24"
                    fill="none"
                    stroke="currentColor"
                    strokeWidth="2"
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  >
                    <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4" />
                    <polyline points="7 10 12 15 17 10" />
                    <line x1="12" y1="15" x2="12" y2="3" />
                  </svg>
                  Download
                </button>
              </div>
              <div className="p-4 flex justify-center">
                <img
                  src={processedImage}
                  alt="Processed Medical Image"
                  className="max-w-full max-h-[70vh] object-contain rounded-lg shadow-md"
                />
              </div>
            </div>
          )}

          {/* Original image display */}
          {!loading && !processedImage && uploadedImagePreview && (
            <div className="w-full max-w-4xl bg-white rounded-lg shadow-md">
              <div className="p-4 border-b">
                <h2 className="text-lg font-bold">Original Image</h2>
                <p className="text-sm text-gray-600">
                  Select preprocessing methods from the sidebar to enhance this image
                </p>
              </div>
              <div className="p-4 flex justify-center">
                <img
                  src={uploadedImagePreview}
                  alt="Original Medical Image"
                  className="max-w-full max-h-[70vh] object-contain rounded-lg shadow-md"
                />
              </div>
            </div>
          )}

          {/* No image uploaded state */}
          {!loading && !uploadedImagePreview && (
            <div className="w-full max-w-md p-6 bg-white rounded-lg shadow-md">
              <div className="flex flex-col items-center text-center">
                <div className="p-6 bg-blue-100 rounded-full mb-4">
                  <svg
                    xmlns="http://www.w3.org/2000/svg"
                    width="48"
                    height="48"
                    viewBox="0 0 24 24"
                    fill="none"
                    stroke="currentColor"
                    strokeWidth="2"
                    strokeLinecap="round"
                    strokeLinejoin="round"
                    className="text-blue-800"
                  >
                    <rect width="18" height="18" x="3" y="3" rx="2" ry="2"></rect>
                    <circle cx="9" cy="9" r="2"></circle>
                    <path d="m21 15-3.086-3.086a2 2 0 0 0-2.828 0L6 21"></path>
                  </svg>
                </div>
                <h3 className="text-xl font-bold mb-2">No Image Uploaded</h3>
                <p className="text-gray-600 mb-4">
                  Please upload a medical image to begin enhancement processing
                </p>
                <p className="text-sm text-gray-500 max-w-xs">
                  Supported formats include DICOM, JPEG, PNG, and TIFF files from various imaging modalities
                </p>
              </div>
            </div>
          )}

          {/* Error state */}
          {error && (
            <div className="w-full max-w-md mt-4 bg-red-50 rounded-lg shadow-md p-4">
              <div className="flex items-center">
                <svg
                  xmlns="http://www.w3.org/2000/svg"
                  width="24"
                  height="24"
                  viewBox="0 0 24 24"
                  fill="none"
                  stroke="currentColor"
                  strokeWidth="2"
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  className="text-red-600 mr-2 flex-shrink-0"
                >
                  <circle cx="12" cy="12" r="10"></circle>
                  <line x1="12" x2="12" y1="8" y2="12"></line>
                  <line x1="12" x2="12.01" y1="16" y2="16"></line>
                </svg>
                <p className="text-red-700">
                  {error}
                </p>
              </div>
            </div>
          )}
        </main>
      </div>
    </div>
  );