'use client';

import React, { useState, useRef, useEffect, useCallback, useMemo } from 'react';
import {
  Input,
  Progress,
  Slider,
  Checkbox,
  Select,
  SelectItem,
  Spacer,
} from '@nextui-org/react';
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
    websocketRef.current = new WebSocket('ws://localhost:8000/imaging/image-processing/ws');

    websocketRef.current.onopen = () => console.log('WebSocket connected.');

    websocketRef.current.onmessage = (event) => {
      if (typeof event.data === 'string') {
        try {
          const data = JSON.parse(event.data);
          if (data.error) {
            setError(data.error);
          }
        } catch (e) {
          setProcessedImage(`data:image/png;base64,${event.data}`);
        }
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
          }, {} as Record<string, boolean>),
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
    debouncedEnhancement(); // Trigger enhancement process
  };

  // Handle active methods change
  const handleMethodChange = (method: string, isChecked: boolean) => {
    setActiveMethods((prev) => {
      const newMethods = isChecked ? [...prev, method] : prev.filter((m) => m !== method);
      debouncedEnhancement();
      return newMethods;
    });
  };

  return (
    <div className="flex flex-col h-screen">
      {/* Main Content */}
      <div className="flex flex-grow">
        {/* Sidebar for Controls */}
        <aside className="w-1/4 p-4 border-r overflow-y-auto">
          <Input
            type="file"
            accept="image/*"
            onChange={handleImageUpload}
            label="Upload Image"
            fullWidth
            css={{ marginBottom: '20px' }}
          />

          <h2 style={{ marginBottom: '10px' }}>Preprocessing Methods</h2>

          {/* Denoising */}
          <div>
            <Checkbox
              isSelected={activeMethods.includes('denoising')}
              onChange={(e) => handleMethodChange('denoising', e.target.checked)}
            >
              Denoising
            </Checkbox>
            {activeMethods.includes('denoising') && (
              <div style={{ paddingLeft: '20px', marginTop: '10px' }}>
                <Select
                  label="Denoise Method"
                  placeholder="Select a method"
                  value={params.denoiseMethod}
                  onChange={(value) => handleParamChange('denoiseMethod', value)}
                  css={{ marginBottom: '10px' }}
                >
                  <SelectItem key="bilateral" value="bilateral">
                    Bilateral
                  </SelectItem>
                  <SelectItem key="gaussian" value="gaussian">
                    Gaussian
                  </SelectItem>
                  <SelectItem key="median" value="median">
                    Median
                  </SelectItem>
                </Select>
                {params.denoiseMethod === 'bilateral' && (
                  <>
                    <Slider
                      label="Diameter"
                      value={[params.bilateralD]}
                      min={1}
                      max={50}
                      step={1}
                      onChange={(value) => handleParamChange('bilateralD', value[0])}
                      css={{ marginBottom: '10px' }}
                    />
                    <Slider
                      label="Sigma Color"
                      value={[params.bilateralSigmaColor]}
                      min={10}
                      max={150}
                      step={1}
                      onChange={(value) => handleParamChange('bilateralSigmaColor', value[0])}
                      css={{ marginBottom: '10px' }}
                    />
                    <Slider
                      label="Sigma Space"
                      value={[params.bilateralSigmaSpace]}
                      min={10}
                      max={150}
                      step={1}
                      onChange={(value) => handleParamChange('bilateralSigmaSpace', value[0])}
                      css={{ marginBottom: '10px' }}
                    />
                  </>
                )}
                {(params.denoiseMethod === 'gaussian' || params.denoiseMethod === 'median') && (
                  <Slider
                    label="Kernel Size"
                    value={[params.gaussianKernelSize]}
                    min={1}
                    max={15}
                    step={2}
                    onChange={(value) => handleParamChange('gaussianKernelSize', value[0])}
                    css={{ marginBottom: '10px' }}
                  />
                )}
                <Spacer y={1} />
              </div>
            )}
          </div>

          {/* Normalization */}
          <div>
            <Checkbox
              isSelected={activeMethods.includes('normalization')}
              onChange={(e) => handleMethodChange('normalization', e.target.checked)}
            >
              Normalization
            </Checkbox>
            {activeMethods.includes('normalization') && (
              <div style={{ paddingLeft: '20px', marginTop: '10px' }}>
                <Slider
                  label="Normalize Min"
                  value={[params.normalizeMin]}
                  min={0}
                  max={255}
                  step={1}
                  onChange={(value) => handleParamChange('normalizeMin', value[0])}
                  css={{ marginBottom: '10px' }}
                />
                <Slider
                  label="Normalize Max"
                  value={[params.normalizeMax]}
                  min={0}
                  max={255}
                  step={1}
                  onChange={(value) => handleParamChange('normalizeMax', value[0])}
                  css={{ marginBottom: '10px' }}
                />
                <Spacer y={1} />
              </div>
            )}
          </div>

          {/* Histogram Equalization */}
          <div>
            <Checkbox
              isSelected={activeMethods.includes('histogramEqualization')}
              onChange={(e) => handleMethodChange('histogramEqualization', e.target.checked)}
            >
              Histogram Equalization
            </Checkbox>
            {/* No adjustable parameters for Histogram Equalization */}
            {activeMethods.includes('histogramEqualization') && (
              <div style={{ paddingLeft: '20px', marginTop: '10px' }}>
                <p>No adjustable parameters.</p>
                <Spacer y={1} />
              </div>
            )}
          </div>

          {/* CLAHE */}
          <div>
            <Checkbox
              isSelected={activeMethods.includes('clahe')}
              onChange={(e) => handleMethodChange('clahe', e.target.checked)}
            >
              CLAHE
            </Checkbox>
            {activeMethods.includes('clahe') && (
              <div style={{ paddingLeft: '20px', marginTop: '10px' }}>
                <Slider
                  label="Clip Limit"
                  value={[params.clipLimit]}
                  min={0.1}
                  max={5.0}
                  step={0.1}
                  onChange={(value) => handleParamChange('clipLimit', value[0])}
                  css={{ marginBottom: '10px' }}
                />
                <Slider
                  label="Tile Grid Size"
                  value={[params.tileGridSize]}
                  min={1}
                  max={16}
                  step={1}
                  onChange={(value) => handleParamChange('tileGridSize', value[0])}
                  css={{ marginBottom: '10px' }}
                />
                <Spacer y={1} />
              </div>
            )}
          </div>

          {/* Unsharp Masking */}
          <div>
            <Checkbox
              isSelected={activeMethods.includes('unsharpMasking')}
              onChange={(e) => handleMethodChange('unsharpMasking', e.target.checked)}
            >
              Unsharp Masking
            </Checkbox>
            {activeMethods.includes('unsharpMasking') && (
              <div style={{ paddingLeft: '20px', marginTop: '10px' }}>
                <Slider
                  label="Unsharp Strength"
                  value={[params.unsharpStrength]}
                  min={0.1}
                  max={3.0}
                  step={0.1}
                  onChange={(value) => handleParamChange('unsharpStrength', value[0])}
                  css={{ marginBottom: '10px' }}
                />
                <Slider
                  label="Gaussian Kernel Size"
                  value={[params.gaussianKernelSize]}
                  min={1}
                  max={15}
                  step={2}
                  onChange={(value) => handleParamChange('gaussianKernelSize', value[0])}
                  css={{ marginBottom: '10px' }}
                />
                <Spacer y={1} />
              </div>
            )}
          </div>

          {/* Windowing */}
          <div>
            <Checkbox
              isSelected={activeMethods.includes('windowing')}
              onChange={(e) => handleMethodChange('windowing', e.target.checked)}
            >
              Windowing
            </Checkbox>
            {activeMethods.includes('windowing') && (
              <div style={{ paddingLeft: '20px', marginTop: '10px' }}>
                <Slider
                  label="Window Width"
                  value={[params.windowWidth]}
                  min={1}
                  max={512}
                  step={1}
                  onChange={(value) => handleParamChange('windowWidth', value[0])}
                  css={{ marginBottom: '10px' }}
                />
                <Slider
                  label="Window Center"
                  value={[params.windowCenter]}
                  min={1}
                  max={512}
                  step={1}
                  onChange={(value) => handleParamChange('windowCenter', value[0])}
                  css={{ marginBottom: '10px' }}
                />
                <Spacer y={1} />
              </div>
            )}
          </div>

          {/* Threshold Segmentation */}
          <div>
            <Checkbox
              isSelected={activeMethods.includes('thresholdSegmentation')}
              onChange={(e) => handleMethodChange('thresholdSegmentation', e.target.checked)}
            >
              Threshold Segmentation
            </Checkbox>
            {activeMethods.includes('thresholdSegmentation') && (
              <div style={{ paddingLeft: '20px', marginTop: '10px' }}>
                <Slider
                  label="Threshold Min"
                  value={[params.thresholdMin]}
                  min={0}
                  max={255}
                  step={1}
                  onChange={(value) => handleParamChange('thresholdMin', value[0])}
                  css={{ marginBottom: '10px' }}
                />
                <Slider
                  label="Threshold Max"
                  value={[params.thresholdMax]}
                  min={0}
                  max={255}
                  step={1}
                  onChange={(value) => handleParamChange('thresholdMax', value[0])}
                  css={{ marginBottom: '10px' }}
                />
                <Spacer y={1} />
              </div>
            )}
          </div>

          {/* Pseudocolor */}
          <div>
            <Checkbox
              isSelected={activeMethods.includes('pseudocolor')}
              onChange={(e) => handleMethodChange('pseudocolor', e.target.checked)}
            >
              Pseudocolor
            </Checkbox>
            {activeMethods.includes('pseudocolor') && (
              <div style={{ paddingLeft: '20px', marginTop: '10px' }}>
                <Select
                  label="Pseudocolor Map"
                  placeholder="Select a map"
                  value={params.pseudocolorMap}
                  onChange={(value) => handleParamChange('pseudocolorMap', value)}
                  css={{ marginBottom: '10px' }}
                >
                  <SelectItem key="COLORMAP_AUTUMN" value="COLORMAP_AUTUMN">
                    Autumn
                  </SelectItem>
                  <SelectItem key="COLORMAP_BONE" value="COLORMAP_BONE">
                    Bone
                  </SelectItem>
                  <SelectItem key="COLORMAP_JET" value="COLORMAP_JET">
                    Jet
                  </SelectItem>
                  <SelectItem key="COLORMAP_WINTER" value="COLORMAP_WINTER">
                    Winter
                  </SelectItem>
                  <SelectItem key="COLORMAP_RAINBOW" value="COLORMAP_RAINBOW">
                    Rainbow
                  </SelectItem>
                  <SelectItem key="COLORMAP_OCEAN" value="COLORMAP_OCEAN">
                    Ocean
                  </SelectItem>
                  <SelectItem key="COLORMAP_SUMMER" value="COLORMAP_SUMMER">
                    Summer
                  </SelectItem>
                  <SelectItem key="COLORMAP_SPRING" value="COLORMAP_SPRING">
                    Spring
                  </SelectItem>
                  <SelectItem key="COLORMAP_COOL" value="COLORMAP_COOL">
                    Cool
                  </SelectItem>
                  <SelectItem key="COLORMAP_HSV" value="COLORMAP_HSV">
                    HSV
                  </SelectItem>
                  <SelectItem key="COLORMAP_PINK" value="COLORMAP_PINK">
                    Pink
                  </SelectItem>
                  <SelectItem key="COLORMAP_HOT" value="COLORMAP_HOT">
                    Hot
                  </SelectItem>
                </Select>
                <Spacer y={1} />
              </div>
            )}
          </div>

          {/* Morphological Operations */}
          <div>
            <Checkbox
              isSelected={activeMethods.includes('morphologicalOperations')}
              onChange={(e) => handleMethodChange('morphologicalOperations', e.target.checked)}
            >
              Morphological Operations
            </Checkbox>
            {activeMethods.includes('morphologicalOperations') && (
              <div style={{ paddingLeft: '20px', marginTop: '10px' }}>
                <Select
                  label="Operation"
                  placeholder="Select an operation"
                  value={params.morphOperation}
                  onChange={(value) => handleParamChange('morphOperation', value)}
                  css={{ marginBottom: '10px' }}
                >
                  <SelectItem key="dilation" value="dilation">
                    Dilation
                  </SelectItem>
                  <SelectItem key="erosion" value="erosion">
                    Erosion
                  </SelectItem>
                  <SelectItem key="opening" value="opening">
                    Opening
                  </SelectItem>
                  <SelectItem key="closing" value="closing">
                    Closing
                  </SelectItem>
                </Select>
                <Slider
                  label="Kernel Size"
                  value={[params.kernelSize]}
                  min={1}
                  max={15}
                  step={2}
                  onChange={(value) => handleParamChange('kernelSize', value[0])}
                  css={{ marginBottom: '10px' }}
                />
                <Spacer y={1} />
              </div>
            )}
          </div>
        </aside>

        {/* Main Image Area */}
        <main className="flex-1 p-4 flex flex-col items-center justify-center">
          {loading && (
            <div className="flex flex-col items-center">
              <Progress indeterminate />
              <p style={{ marginTop: '10px' }}>Processing image...</p>
            </div>
          )}
          {!loading && processedImage && (
            <img
              src={processedImage}
              alt="Processed"
              style={{
                maxWidth: '100%',
                borderRadius: '10px',
                boxShadow: '0 4px 6px rgba(0, 0, 0, 0.1)',
              }}
            />
          )}
          {!loading && !processedImage && uploadedImagePreview && (
            <img
              src={uploadedImagePreview}
              alt="Uploaded"
              style={{
                maxWidth: '100%',
                borderRadius: '10px',
                boxShadow: '0 4px 6px rgba(0, 0, 0, 0.1)',
              }}
            />
          )}
          {error && (
            <p style={{ color: 'red', marginTop: '20px' }}>
              {error}
            </p>
          )}
        </main>
      </div>
    </div>
  );
}
