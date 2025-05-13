'use client';

import React, { useState, useRef, useEffect, useCallback } from 'react';
import {
  Input,
  Progress,
  Slider,
  Checkbox,
  Select,
  SelectItem,
  Spacer,
  Card,
  CardBody,
  CardHeader,
  Divider,
  Button,
  Tooltip,
} from '@nextui-org/react';

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
    bilateralD: 9, // Diameter of each pixel neighborhood
    bilateralSigmaColor: 75, // Filter sigma in the color space
    bilateralSigmaSpace: 75, // Filter sigma in the coordinate space
    gaussianKernelSize: 5, // Gaussian kernel size (must be odd)
    medianKernelSize: 5, // Median filter kernel size (must be odd)

    // CLAHE parameters
    clipLimit: 2.0, // Threshold for contrast limiting
    tileGridSize: 8, // Size of grid for histogram equalization

    // Sharpening parameters
    unsharpStrength: 1.5, // Strength of unsharp mask

    // Intensity transformation parameters
    normalizeMin: 0, // Lower bound for normalization
    normalizeMax: 255, // Upper bound for normalization

    // Windowing parameters (for radiological images)
    windowWidth: 255, // Window width for contrast adjustment
    windowCenter: 127, // Window center for brightness adjustment

    // Segmentation parameters
    thresholdMin: 50, // Lower threshold boundary
    thresholdMax: 200, // Upper threshold boundary

    // Visualization parameters
    pseudocolorMap: 'COLORMAP_JET', // Color mapping for pseudocolor

    // Morphological parameters
    kernelSize: 3, // Structuring element size (must be odd)
    morphOperation: 'dilation', // Type of morphological operation
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
    // Establish connection to backend service
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
            // If not JSON, expect base64 image data
            setProcessedImage(`data:image/png;base64,${event.data}`);
            setLoading(false);
          }
        }
      };

      websocketRef.current.onerror = (err) => {
        // setError('WebSocket communication error. Please check server connection.');
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

    // Clean up connection on component unmount
    return () => {
      if (websocketRef.current) {
        websocketRef.current.close();
      }
    };
  }, []);

  /**
   * Initiates the image enhancement process by sending image data and parameters
   * to the backend processing service via WebSocket.
   */
  const startEnhancement = useCallback(() => {
    if (!image) {
      setError('Please upload a medical image to process.');
      return;
    }

    if (!websocketRef.current || websocketRef.current.readyState !== WebSocket.OPEN) {
      setError('WebSocket connection is not established. Please refresh the page.');
      return;
    }

    // If no active methods, show original image
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

        // Send image data and processing parameters to server
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

  // Apply processing when parameters or methods change
  useEffect(() => {
    if (image) {
      // Add debounce to avoid excessive processing
      const debounceTimeout = setTimeout(() => {
        startEnhancement();
      }, 500);

      return () => clearTimeout(debounceTimeout);
    }
  }, [image, params, activeMethods, startEnhancement]);

  /**
   * Handles image file selection and validates the uploaded file
   * @param {React.ChangeEvent<HTMLInputElement>} event - The file input change event
   */
  const handleImageUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const selectedFile = event.target.files?.[0];

    if (selectedFile) {
      // Validate file type
      const validImageTypes = ['image/jpeg', 'image/png', 'image/dicom', 'image/tiff'];
      if (!validImageTypes.includes(selectedFile.type) && !selectedFile.name.endsWith('.dcm')) {
        setError('Please upload a valid medical image format (JPEG, PNG, DICOM, TIFF).');
        return;
      }

      // Create preview and store image
      const objectUrl = URL.createObjectURL(selectedFile);
      setImage(selectedFile);
      setUploadedImagePreview(objectUrl);
      setOriginalImage(objectUrl);
      setProcessedImage(null);
      setError(null);

      // Reset active methods when new image is uploaded
      setActiveMethods([]);
    }
  };

  /**
   * Updates a specific parameter value in the parameters state
   * @param {string} name - The parameter name to update
   * @param {any} value - The new parameter value
   */
  const handleParamChange = (name: string, value: any) => {
    setParams((prev) => ({
      ...prev,
      [name]: value,
    }));
  };

  /**
   * Toggles an image processing method on or off
   * @param {string} method - The method name to toggle
   * @param {boolean} isChecked - Whether the method should be active
   */
  const handleMethodChange = (method: string, isChecked: boolean) => {
    setActiveMethods((prev) => {
      return isChecked
        ? [...prev, method]
        : prev.filter((m) => m !== method);
    });
  };

  /**
   * Resets all parameters to their default values
   */
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

  /**
   * Toggles all processing methods on or off
   * @param {boolean} enable - Whether to enable or disable all methods
   */
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

      {/* Main Content */}
      <div className="flex flex-grow overflow-hidden">
        {/* Sidebar for Controls */}
        <aside className="w-1/3 p-4 border-r overflow-y-auto bg-white">
          <Card className="mb-4">
            <CardHeader className="pb-0 pt-2 px-4 flex-col items-start">
              <h2 className="text-lg font-bold">Image Input</h2>
            </CardHeader>
            <CardBody>
              <Input
                type="file"
                accept="image/*,.dcm"
                onChange={handleImageUpload}
                label="Upload Medical Image"
                description="Supported formats: JPEG, PNG, DICOM, TIFF"
                fullWidth
              />

              <div className="flex space-x-2 mt-4">
                <Button
                  size="sm"
                  color="primary"
                  onClick={() => toggleAllMethods(true)}
                  disabled={!image}
                >
                  Enable All Methods
                </Button>
                <Button
                  size="sm"
                  color="danger"
                  onClick={() => toggleAllMethods(false)}
                  disabled={!image}
                >
                  Disable All Methods
                </Button>
                <Button
                  size="sm"
                  color="warning"
                  onClick={handleResetParams}
                  disabled={!image}
                >
                  Reset Parameters
                </Button>
              </div>
            </CardBody>
          </Card>

          <Card className="mb-4">
            <CardHeader className="pb-0 pt-2 px-4 flex-col items-start">
              <h2 className="text-lg font-bold">Processing Methods</h2>
              <p className="text-sm text-gray-500">Select techniques to enhance the image</p>
            </CardHeader>
            <CardBody>
              {/* Denoising */}
              <div className="mb-4">
                <div className="flex items-center">
                  <Checkbox
                    isSelected={activeMethods.includes('denoising')}
                    onChange={(e) => handleMethodChange('denoising', e.target.checked)}
                  >
                    <span className="font-medium">Denoising</span>
                  </Checkbox>
                  <Tooltip content={methodDescriptions.denoising}>
                    <span className="ml-2 text-gray-500 cursor-help text-sm">ⓘ</span>
                  </Tooltip>
                </div>

                {activeMethods.includes('denoising') && (
                  <div className="pl-6 mt-2">
                    <Select
                      label="Filter Algorithm"
                      placeholder="Select algorithm"
                      value={params.denoiseMethod}
                      onChange={(e) => handleParamChange('denoiseMethod', e.target.value)}
                      className="mb-2"
                    >
                      <SelectItem key="bilateral" value="bilateral">
                        Bilateral Filter
                      </SelectItem>
                      <SelectItem key="gaussian" value="gaussian">
                        Gaussian Filter
                      </SelectItem>
                      <SelectItem key="median" value="median">
                        Median Filter
                      </SelectItem>
                    </Select>

                    {params.denoiseMethod === 'bilateral' && (
                      <>
                        <Slider
                          label={`Diameter: ${params.bilateralD}`}
                          value={[params.bilateralD]}
                          min={3}
                          max={25}
                          step={2}
                          onChange={(value) => handleParamChange('bilateralD', value[0])}
                          className="mb-2"
                        />
                        <Slider
                          label={`Color Sigma: ${params.bilateralSigmaColor}`}
                          value={[params.bilateralSigmaColor]}
                          min={10}
                          max={150}
                          step={5}
                          onChange={(value) => handleParamChange('bilateralSigmaColor', value[0])}
                          className="mb-2"
                        />
                        <Slider
                          label={`Space Sigma: ${params.bilateralSigmaSpace}`}
                          value={[params.bilateralSigmaSpace]}
                          min={10}
                          max={150}
                          step={5}
                          onChange={(value) => handleParamChange('bilateralSigmaSpace', value[0])}
                          className="mb-2"
                        />
                      </>
                    )}

                    {params.denoiseMethod === 'gaussian' && (
                      <Slider
                        label={`Kernel Size: ${params.gaussianKernelSize}`}
                        value={[params.gaussianKernelSize]}
                        min={1}
                        max={15}
                        step={2}
                        onChange={(value) => handleParamChange('gaussianKernelSize', value[0])}
                        className="mb-2"
                      />
                    )}

                    {params.denoiseMethod === 'median' && (
                      <Slider
                        label={`Kernel Size: ${params.medianKernelSize}`}
                        value={[params.medianKernelSize]}
                        min={1}
                        max={15}
                        step={2}
                        onChange={(value) => handleParamChange('medianKernelSize', value[0])}
                        className="mb-2"
                      />
                    )}
                  </div>
                )}
              </div>

              <Divider className="my-2" />

              {/* Normalization */}
              <div className="mb-4">
                <div className="flex items-center">
                  <Checkbox
                    isSelected={activeMethods.includes('normalization')}
                    onChange={(e) => handleMethodChange('normalization', e.target.checked)}
                  >
                    <span className="font-medium">Intensity Normalization</span>
                  </Checkbox>
                  <Tooltip content={methodDescriptions.normalization}>
                    <span className="ml-2 text-gray-500 cursor-help text-sm">ⓘ</span>
                  </Tooltip>
                </div>

                {activeMethods.includes('normalization') && (
                  <div className="pl-6 mt-2">
                    <Slider
                      label={`Min Intensity: ${params.normalizeMin}`}
                      value={[params.normalizeMin]}
                      min={0}
                      max={255}
                      step={1}
                      onChange={(value) => handleParamChange('normalizeMin', value[0])}
                      className="mb-2"
                    />
                    <Slider
                      label={`Max Intensity: ${params.normalizeMax}`}
                      value={[params.normalizeMax]}
                      min={0}
                      max={255}
                      step={1}
                      onChange={(value) => handleParamChange('normalizeMax', value[0])}
                      className="mb-2"
                    />
                  </div>
                )}
              </div>

              <Divider className="my-2" />

              {/* Histogram Equalization */}
              <div className="mb-4">
                <div className="flex items-center">
                  <Checkbox
                    isSelected={activeMethods.includes('histogramEqualization')}
                    onChange={(e) => handleMethodChange('histogramEqualization', e.target.checked)}
                  >
                    <span className="font-medium">Histogram Equalization</span>
                  </Checkbox>
                  <Tooltip content={methodDescriptions.histogramEqualization}>
                    <span className="ml-2 text-gray-500 cursor-help text-sm">ⓘ</span>
                  </Tooltip>
                </div>

                {activeMethods.includes('histogramEqualization') && (
                  <div className="pl-6 mt-2">
                    <p className="text-sm text-gray-600">Global histogram equalization with no adjustable parameters.</p>
                  </div>
                )}
              </div>

              <Divider className="my-2" />

              {/* CLAHE */}
              <div className="mb-4">
                <div className="flex items-center">
                  <Checkbox
                    isSelected={activeMethods.includes('clahe')}
                    onChange={(e) => handleMethodChange('clahe', e.target.checked)}
                  >
                    <span className="font-medium">CLAHE</span>
                  </Checkbox>
                  <Tooltip content={methodDescriptions.clahe}>
                    <span className="ml-2 text-gray-500 cursor-help text-sm">ⓘ</span>
                  </Tooltip>
                </div>

                {activeMethods.includes('clahe') && (
                  <div className="pl-6 mt-2">
                    <Slider
                      label={`Clip Limit: ${params.clipLimit.toFixed(1)}`}
                      value={[params.clipLimit]}
                      min={0.5}
                      max={8.0}
                      step={0.1}
                      onChange={(value) => handleParamChange('clipLimit', value[0])}
                      className="mb-2"
                    />
                    <Slider
                      label={`Tile Grid Size: ${params.tileGridSize}`}
                      value={[params.tileGridSize]}
                      min={2}
                      max={16}
                      step={2}
                      onChange={(value) => handleParamChange('tileGridSize', value[0])}
                      className="mb-2"
                    />
                  </div>
                )}
              </div>

              <Divider className="my-2" />

              {/* Unsharp Masking */}
              <div className="mb-4">
                <div className="flex items-center">
                  <Checkbox
                    isSelected={activeMethods.includes('unsharpMasking')}
                    onChange={(e) => handleMethodChange('unsharpMasking', e.target.checked)}
                  >
                    <span className="font-medium">Edge Enhancement</span>
                  </Checkbox>
                  <Tooltip content={methodDescriptions.unsharpMasking}>
                    <span className="ml-2 text-gray-500 cursor-help text-sm">ⓘ</span>
                  </Tooltip>
                </div>

                {activeMethods.includes('unsharpMasking') && (
                  <div className="pl-6 mt-2">
                    <Slider
                      label={`Strength: ${params.unsharpStrength.toFixed(1)}`}
                      value={[params.unsharpStrength]}
                      min={0.1}
                      max={5.0}
                      step={0.1}
                      onChange={(value) => handleParamChange('unsharpStrength', value[0])}
                      className="mb-2"
                    />
                    <Slider
                      label={`Blur Radius: ${params.gaussianKernelSize}`}
                      value={[params.gaussianKernelSize]}
                      min={3}
                      max={15}
                      step={2}
                      onChange={(value) => handleParamChange('gaussianKernelSize', value[0])}
                      className="mb-2"
                    />
                  </div>
                )}
              </div>

              <Divider className="my-2" />

              {/* Windowing */}
              <div className="mb-4">
                <div className="flex items-center">
                  <Checkbox
                    isSelected={activeMethods.includes('windowing')}
                    onChange={(e) => handleMethodChange('windowing', e.target.checked)}
                  >
                    <span className="font-medium">Radiological Windowing</span>
                  </Checkbox>
                  <Tooltip content={methodDescriptions.windowing}>
                    <span className="ml-2 text-gray-500 cursor-help text-sm">ⓘ</span>
                  </Tooltip>
                </div>

                {activeMethods.includes('windowing') && (
                  <div className="pl-6 mt-2">
                    <Slider
                      label={`Window Width: ${params.windowWidth}`}
                      value={[params.windowWidth]}
                      min={1}
                      max={512}
                      step={1}
                      onChange={(value) => handleParamChange('windowWidth', value[0])}
                      className="mb-2"
                    />
                    <Slider
                      label={`Window Center: ${params.windowCenter}`}
                      value={[params.windowCenter]}
                      min={1}
                      max={512}
                      step={1}
                      onChange={(value) => handleParamChange('windowCenter', value[0])}
                      className="mb-2"
                    />
                  </div>
                )}
              </div>

              <Divider className="my-2" />

              {/* Threshold Segmentation */}
              <div className="mb-4">
                <div className="flex items-center">
                  <Checkbox
                    isSelected={activeMethods.includes('thresholdSegmentation')}
                    onChange={(e) => handleMethodChange('thresholdSegmentation', e.target.checked)}
                  >
                    <span className="font-medium">Threshold Segmentation</span>
                  </Checkbox>
                  <Tooltip content={methodDescriptions.thresholdSegmentation}>
                    <span className="ml-2 text-gray-500 cursor-help text-sm">ⓘ</span>
                  </Tooltip>
                </div>

                {activeMethods.includes('thresholdSegmentation') && (
                  <div className="pl-6 mt-2">
                    <Slider
                      label={`Lower Threshold: ${params.thresholdMin}`}
                      value={[params.thresholdMin]}
                      min={0}
                      max={255}
                      step={1}
                      onChange={(value) => handleParamChange('thresholdMin', value[0])}
                      className="mb-2"
                    />
                    <Slider
                      label={`Upper Threshold: ${params.thresholdMax}`}
                      value={[params.thresholdMax]}
                      min={0}
                      max={255}
                      step={1}
                      onChange={(value) => handleParamChange('thresholdMax', value[0])}
                      className="mb-2"
                    />
                  </div>
                )}
              </div>

              <Divider className="my-2" />

              {/* Pseudocolor */}
              <div className="mb-4">
                <div className="flex items-center">
                  <Checkbox
                    isSelected={activeMethods.includes('pseudocolor')}
                    onChange={(e) => handleMethodChange('pseudocolor', e.target.checked)}
                  >
                    <span className="font-medium">Pseudocolor Mapping</span>
                  </Checkbox>
                  <Tooltip content={methodDescriptions.pseudocolor}>
                    <span className="ml-2 text-gray-500 cursor-help text-sm">ⓘ</span>
                  </Tooltip>
                </div>

                {activeMethods.includes('pseudocolor') && (
                  <div className="pl-6 mt-2">
                    <Select
                      label="Color Map"
                      placeholder="Select color map"
                      value={params.pseudocolorMap}
                      onChange={(e) => handleParamChange('pseudocolorMap', e.target.value)}
                      className="mb-2"
                    >
                      <SelectItem key="COLORMAP_JET" value="COLORMAP_JET">Jet</SelectItem>
                      <SelectItem key="COLORMAP_VIRIDIS" value="COLORMAP_VIRIDIS">Viridis</SelectItem>
                      <SelectItem key="COLORMAP_HOT" value="COLORMAP_HOT">Hot</SelectItem>
                      <SelectItem key="COLORMAP_BONE" value="COLORMAP_BONE">Bone</SelectItem>
                      <SelectItem key="COLORMAP_PLASMA" value="COLORMAP_PLASMA">Plasma</SelectItem>
                      <SelectItem key="COLORMAP_INFERNO" value="COLORMAP_INFERNO">Inferno</SelectItem>
                      <SelectItem key="COLORMAP_MAGMA" value="COLORMAP_MAGMA">Magma</SelectItem>
                      <SelectItem key="COLORMAP_CIVIDIS" value="COLORMAP_CIVIDIS">Cividis</SelectItem>
                      <SelectItem key="COLORMAP_RAINBOW" value="COLORMAP_RAINBOW">Rainbow</SelectItem>
                      <SelectItem key="COLORMAP_OCEAN" value="COLORMAP_OCEAN">Ocean</SelectItem>
                      <SelectItem key="COLORMAP_TURBO" value="COLORMAP_TURBO">Turbo</SelectItem>
                    </Select>
                  </div>
                )}
              </div>

              <Divider className="my-2" />

              {/* Morphological Operations */}
              <div className="mb-4">
                <div className="flex items-center">
                  <Checkbox
                    isSelected={activeMethods.includes('morphologicalOperations')}
                    onChange={(e) => handleMethodChange('morphologicalOperations', e.target.checked)}
                  >
                    <span className="font-medium">Morphological Operations</span>
                  </Checkbox>
                  <Tooltip content={methodDescriptions.morphologicalOperations}>
                    <span className="ml-2 text-gray-500 cursor-help text-sm">ⓘ</span>
                  </Tooltip>
                </div>

                {activeMethods.includes('morphologicalOperations') && (
                  <div className="pl-6 mt-2">
                    <Select
                      label="Operation Type"
                      placeholder="Select operation"
                      value={params.morphOperation}
                      onChange={(e) => handleParamChange('morphOperation', e.target.value)}
                      className="mb-2"
                    >
                      <SelectItem key="dilation" value="dilation">Dilation</SelectItem>
                      <SelectItem key="erosion" value="erosion">Erosion</SelectItem>
                      <SelectItem key="opening" value="opening">Opening</SelectItem>
                      <SelectItem key="closing" value="closing">Closing</SelectItem>
                      <SelectItem key="gradient" value="gradient">Morphological Gradient</SelectItem>
                      <SelectItem key="tophat" value="tophat">Top Hat</SelectItem>
                      <SelectItem key="blackhat" value="blackhat">Black Hat</SelectItem>
                    </Select>
                    <Slider
                      label={`Kernel Size: ${params.kernelSize}`}
                      value={[params.kernelSize]}
                      min={1}
                      max={15}
                      step={2}
                      onChange={(value) => handleParamChange('kernelSize', value[0])}
                      className="mb-2"
                    />
                  </div>
                )}
              </div>
            </CardBody>
          </Card>
        </aside>

        {/* Main Image Display Area */}
        <main className="flex-1 p-4 flex flex-col items-center justify-center bg-gray-100">
          {/* Processing indicator */}
          {loading && (
            <Card className="w-full max-w-md p-4">
              <CardBody className="flex flex-col items-center">
                <Progress
                  indeterminated
                  color="primary"
                  className="w-full mb-4"
                  aria-label="Processing image..."
                />
                <p className="text-center text-gray-700">
                  Applying image processing algorithms...
                </p>
              </CardBody>
            </Card>
          )}

          {/* Processed image display */}
          {!loading && processedImage && (
            <Card className="w-full max-w-4xl">
              <CardHeader className="pb-0 pt-2 px-4">
                <h2 className="text-lg font-bold">Enhanced Image</h2>
                <p className="text-sm text-gray-600">
                  Applied methods: {activeMethods.length > 0 ?
                    activeMethods.map(method => method.charAt(0).toUpperCase() + method.slice(1)).join(', ') :
                    'None (showing original)'}
                </p>
              </CardHeader>
              <CardBody className="flex justify-center py-4">
                <img
                  src={processedImage}
                  alt="Processed Medical Image"
                  className="max-w-full max-h-[70vh] object-contain rounded-lg shadow-md"
                />
              </CardBody>
            </Card>
          )}

          {/* Original image display */}
          {!loading && !processedImage && uploadedImagePreview && (
            <Card className="w-full max-w-4xl">
              <CardHeader className="pb-0 pt-2 px-4">
                <h2 className="text-lg font-bold">Original Image</h2>
                <p className="text-sm text-gray-600">
                  Select preprocessing methods from the sidebar to enhance this image
                </p>
              </CardHeader>
              <CardBody className="flex justify-center py-4">
                <img
                  src={uploadedImagePreview}
                  alt="Original Medical Image"
                  className="max-w-full max-h-[70vh] object-contain rounded-lg shadow-md"
                />
              </CardBody>
            </Card>
          )}

          {/* No image uploaded state */}
          {!loading && !uploadedImagePreview && (
            <Card className="w-full max-w-md p-6 bg-gray-50">
              <CardBody className="flex flex-col items-center text-center">
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
              </CardBody>
            </Card>
          )}

          {/* Error state */}
          {error && (
            <Card className="w-full max-w-md mt-4 bg-red-50">
              <CardBody className="flex items-center">
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
                  className="text-red-600 mr-2"
                >
                  <circle cx="12" cy="12" r="10"></circle>
                  <line x1="12" x2="12" y1="8" y2="12"></line>
                  <line x1="12" x2="12.01" y1="16" y2="16"></line>
                </svg>
                <p className="text-red-700">
                  {error}
                </p>
              </CardBody>
            </Card>
          )}
        </main>
      </div>

      {/* Footer with academic citations */}
      {/*<footer className="p-4 bg-gray-200 text-sm text-gray-600">*/}
      {/*  <div className="max-w-5xl mx-auto">*/}
      {/*    <h3 className="font-bold mb-1">References:</h3>*/}
      {/*    <ul className="list-disc pl-5">*/}
      {/*      <li>Gonzalez R.C., Woods R.E. (2018). "Digital Image Processing." 4th Edition, Pearson.</li>*/}
      {/*      <li>Pizer S.M., et al. (1987). "Adaptive Histogram Equalization and Its Variations." Computer Vision, Graphics, and Image Processing, 39(3), 355-368.</li>*/}
      {/*      <li>Bankman I.N. (2008). "Handbook of Medical Image Processing and Analysis." 2nd Edition, Academic Press.</li>*/}
      {/*      <li>Wang Z., Bovik A.C. (2009). "Mean squared error: Love it or leave it?" IEEE Signal Processing Magazine, 26(1), 98-117.</li>*/}
      {/*    </ul>*/}
      {/*  </div>*/}
      {/*</footer>*/}
    </div>
  );
}