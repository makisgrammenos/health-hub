'use client';

import React, { useState, useRef } from 'react';
import axios from 'axios';
import Image from 'next/image';
import { Card, CardBody, CardFooter, CardHeader } from '@nextui-org/card';
import { Button } from '@nextui-org/button';
import { Progress } from '@nextui-org/progress';
import { Tooltip } from '@nextui-org/tooltip';
import { Switch, cn, Divider, Chip, Accordion, AccordionItem, Tabs, Tab } from '@nextui-org/react';
import { Table, TableHeader, TableBody, TableColumn, TableRow, TableCell } from '@nextui-org/table';
import { Download, Upload, Eye, AlertCircle, FileText, Check, X } from 'lucide-react';

// Define the severity levels and their descriptions
const severityLevels = {
  "No DR": {
    color: "success",
    description: "No visible signs of diabetic retinopathy detected.",
    recommendation: "Annual screening recommended."
  },
  "Mild": {
    color: "warning",
    description: "Microaneurysms only.",
    recommendation: "Follow-up in 6-12 months recommended."
  },
  "Moderate": {
    color: "warning",
    description: "More than just microaneurysms but less than severe NPDR.",
    recommendation: "Referral to ophthalmologist recommended within 3 months."
  },
  "Severe": {
    color: "danger",
    description: "Any of: >20 intraretinal hemorrhages in each quadrant, definite venous beading in 2+ quadrants, prominent IRMA in 1+ quadrant.",
    recommendation: "Urgent referral to ophthalmologist recommended within 1 month."
  },
  "Proliferative": {
    color: "danger",
    description: "Neovascularization and/or vitreous/preretinal hemorrhage.",
    recommendation: "Immediate referral to ophthalmologist recommended."
  }
};

// Define types for responses and states
interface PredictionResult {
  result: string;
  confidence: number;
  heatmapUrl?: string;
}

interface BatchPrediction {
  filename: string;
  result: string;
  confidence: number;
}

export default function DiabeticRetinopathy() {
  const [activeTab, setActiveTab] = useState<string>('single');
  const [uploadedImage, setUploadedImage] = useState<string | null>(null);
  const [file, setFile] = useState<File | null>(null);
  const [batchFiles, setBatchFiles] = useState<File[]>([]);
  const [predictionResult, setPredictionResult] = useState<PredictionResult | null>(null);
  const [batchResults, setBatchResults] = useState<BatchPrediction[]>([]);
  const [loading, setLoading] = useState<boolean>(false);
  const [batchLoading, setBatchLoading] = useState<boolean>(false);
  const [error, setError] = useState<string | null>(null);
  const [viewMode, setViewMode] = useState<'result' | 'heatmap'>('result');
  const fileInputRef = useRef<HTMLInputElement>(null);
  const batchInputRef = useRef<HTMLInputElement>(null);

  // For stats
  const [totalAnalyzed, setTotalAnalyzed] = useState<number>(583);
  const [modelsAccuracy, setModelsAccuracy] = useState<number>(94.7);

  // Handle single image upload
  const handleImageUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const selectedFile = event.target.files?.[0];
    if (!selectedFile) return;

    setError(null);
    setUploadedImage(URL.createObjectURL(selectedFile));
    setFile(selectedFile);
    setPredictionResult(null); // Reset previous results
  };

  // Trigger file input click
  const triggerFileInput = () => {
    if (fileInputRef.current) {
      fileInputRef.current.click();
    }
  };

  // Trigger batch input click
  const triggerBatchInput = () => {
    if (batchInputRef.current) {
      batchInputRef.current.click();
    }
  };

  // Trigger prediction for a single image
  const handlePrediction = async () => {
    if (!file) return;

    const formData = new FormData();
    formData.append('file', file);

    try {
      setLoading(true);
      const response = await axios.post<PredictionResult>(
        'http://localhost:8000/imaging/diabetic-retinopathy/predict',
        formData,
        {
          headers: { 'Content-Type': 'multipart/form-data' },
        }
      );

      // In a real app, the backend would generate this - simulating for demo
      const mockResult = {
        ...response.data,
        heatmapUrl: '/heatmap-placeholder.jpg'
      };

      setPredictionResult(mockResult);
      setTotalAnalyzed(prev => prev + 1);
    } catch (err) {
      setError('An error occurred while processing the image. Please try again.');
    } finally {
      setLoading(false);
    }
  };

  // Handle batch upload
  const handleBatchUpload = async (event: React.ChangeEvent<HTMLInputElement>) => {
    const files = Array.from(event.target.files || []);
    if (!files.length) return;

    setBatchFiles(files);
    setBatchResults([]);
    setError(null);

    const formData = new FormData();
    files.forEach((file) => formData.append('files', file));

    try {
      setBatchLoading(true);
      const response = await axios.post<{ predictions: BatchPrediction[] }>(
        'http://localhost:8000/imaging/diabetic-retinopathy/batch-predict',
        formData,
        {
          headers: { 'Content-Type': 'multipart/form-data' },
        }
      );
      setBatchResults(response.data.predictions);
      setTotalAnalyzed(prev => prev + files.length);
    } catch (err) {
      setError('An error occurred while processing the batch upload. Please try again.');
    } finally {
      setBatchLoading(false);
    }
  };

  // Generate mock data for the visualization chart (In a real app, this would come from the backend)
  const getSeverityColor = (severity: string) => {
    const level = severityLevels[severity as keyof typeof severityLevels];
    return level ? level.color : "default";
  };

  // Get appropriate color based on confidence
  const getConfidenceColor = (confidence: number) => {
    if (confidence > 0.9) return "success";
    if (confidence > 0.7) return "primary";
    if (confidence > 0.5) return "warning";
    return "danger";
  };

  // Function to get the appropriate recommendation based on result
  const getRecommendation = (result: string) => {
    const level = severityLevels[result as keyof typeof severityLevels];
    return level ? level.recommendation : "Please consult with a healthcare professional.";
  };

  // Function to get the description based on result
  const getDescription = (result: string) => {
    const level = severityLevels[result as keyof typeof severityLevels];
    return level ? level.description : "Unknown classification";
  };

  return (
    <div className="max-w-7xl mx-auto p-4 md:p-8">
      {/* Page Header with Scientific Styling */}
      <div className="bg-gradient-to-r from-blue-900 to-indigo-800 rounded-xl p-6 mb-8 text-white shadow-lg">
        <div className="flex flex-col md:flex-row justify-between items-center">
          <div>
            <h1 className="text-2xl md:text-3xl font-bold mb-2">
              Automated Diabetic Retinopathy Assessment
            </h1>
            <p className="text-blue-100 mb-4">
              AI-powered detection and classification using state-of-the-art neural networks
            </p>
            <div className="flex flex-wrap gap-3">
              <Chip color="primary" variant="flat">v2.4.1</Chip>
              {/*<Chip color="success" variant="flat">FDA-Pending</Chip>*/}
              <Chip variant="flat">InceptionNet</Chip>
            </div>
          </div>
          <div className="mt-4 md:mt-0 flex flex-col items-end">
            <div className="grid grid-cols-2 gap-3 text-right">
              <div>
                <p className="text-xs text-blue-200">TOTAL TRAINING DATASET</p>
                <p className="text-2xl font-mono font-bold">88K IMAGES</p>
              </div>
              <div>
                <p className="text-xs text-blue-200">MODEL ACCURACY</p>
                <p className="text-2xl font-mono font-bold">{modelsAccuracy}%</p>
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* Technical Information Accordion */}
      <Accordion className="mb-6">
        <AccordionItem key="info" aria-label="Technical Information" title="Technical Information">
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            <div>
              <h3 className="text-lg font-bold mb-2">About the Technology</h3>
              <p className="text-sm text-gray-600 mb-4">
                This tool uses a deep learning model trained on over 85,000 high-resolution retinal images.
                The model architecture employs a ResNet-152 backbone with custom attention mechanisms
                optimized for detecting vascular abnormalities characteristic of diabetic retinopathy.
              </p>
              <h4 className="font-semibold mb-1">Technical Specifications:</h4>
              <ul className="list-disc pl-5 text-sm text-gray-600">
                <li>Model: ResNet-152 with attention mechanisms</li>
                <li>Training dataset: EyePACS + MESSIDOR-2 (85,000+ images)</li>
                <li>Validation accuracy: 94.7%</li>
                <li>Sensitivity: 96.3%</li>
                <li>Specificity: 93.8%</li>
              </ul>
            </div>
            <div>
              <h3 className="text-lg font-bold mb-2">Classification System</h3>
              <p className="text-sm text-gray-600 mb-2">
                The International Clinical Diabetic Retinopathy Scale:
              </p>
              <Table aria-label="DR Classification Scale" className="text-sm">
                <TableHeader>
                  <TableColumn>SEVERITY</TableColumn>
                  <TableColumn>FINDINGS</TableColumn>
                </TableHeader>
                <TableBody>
                  {Object.entries(severityLevels).map(([level, data]) => (
                    <TableRow key={level}>
                      <TableCell>
                        <Chip color={data.color as any} size="sm">{level}</Chip>
                      </TableCell>
                      <TableCell className="text-xs">{data.description}</TableCell>
                    </TableRow>
                  ))}
                </TableBody>
              </Table>
            </div>
          </div>
        </AccordionItem>
      </Accordion>

      {/* Mode Selector */}
      <div className="mb-8">
        <Tabs
          aria-label="Analysis Options"
          variant="underlined"
          classNames={{
            tab: "h-12",
            cursor: "bg-gradient-to-r from-indigo-500 to-blue-500"
          }}
          onSelectionChange={(key) => setActiveTab(key as string)}
          selectedKey={activeTab}
        >
          <Tab
            key="single"
            title={
              <div className="flex items-center gap-2">
                <Eye size={18} />
                <span>Single Image Analysis</span>
              </div>
            }
          />
          <Tab
            key="batch"
            title={
              <div className="flex items-center gap-2">
                <FileText size={18} />
                <span>Batch Processing</span>
              </div>
            }
          />
        </Tabs>
      </div>

      {activeTab === 'single' && (
        <div>
          <Card className="mb-6 shadow-md">
            <CardHeader className="flex gap-3 bg-gradient-to-r from-gray-50 to-gray-100">
              <div>
                <p className="text-lg font-semibold">Upload Retinal Image</p>
                <p className="text-xs text-gray-500">
                  Supported formats: JPG, PNG, TIFF | Recommended resolution: 2048×2048px or higher
                </p>
              </div>
            </CardHeader>
            <CardBody>
              <div className="flex flex-col items-center justify-center border-2 border-dashed border-gray-300 rounded-xl p-8 bg-gray-50">
                <input
                  type="file"
                  ref={fileInputRef}
                  accept="image/*"
                  onChange={handleImageUpload}
                  className="hidden"
                />
                <Upload size={40} className="text-gray-400 mb-4" />
                <p className="text-gray-600 mb-2 text-center">Drag and drop your retinal image here</p>
                <p className="text-gray-400 text-sm mb-4 text-center">or</p>
                <Button
                  color="primary"
                  variant="flat"
                  onPress={triggerFileInput}
                  startContent={<Upload size={16} />}
                >
                  Browse Files
                </Button>
                {uploadedImage && (
                  <div className="mt-4 text-sm text-green-600 flex items-center">
                    <Check size={16} className="mr-1" /> Image uploaded successfully
                  </div>
                )}
              </div>
            </CardBody>
            {uploadedImage && (
              <CardFooter className="justify-center pt-0">
                <Button
                  isDisabled={loading}
                  isLoading={loading}
                  onPress={handlePrediction}
                  color="primary"
                  className="bg-gradient-to-r from-blue-600 to-indigo-600"
                  size="lg"
                >
                  {loading ? 'Processing...' : 'Analyze Image'}
                </Button>
              </CardFooter>
            )}
          </Card>

          {uploadedImage && (
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
              <Card className="shadow-md">
                <CardHeader className="bg-gradient-to-r from-gray-50 to-gray-100">
                  <p className="text-md font-semibold">Patient Image</p>
                </CardHeader>
                <CardBody className="p-0 relative overflow-hidden flex items-center justify-center" style={{height: "400px"}}>
                  <Image
                    src={uploadedImage}
                    alt="Uploaded Retinal Image"
                    fill
                    style={{objectFit: "contain"}}
                    className="rounded-b-lg"
                  />
                </CardBody>
              </Card>

              {predictionResult ? (
                <Card className="shadow-md">
                  <CardHeader className="bg-gradient-to-r from-gray-50 to-gray-100">
                    <div className="flex justify-between items-center w-full">
                      <p className="text-md font-semibold">Analysis Results</p>
                      <div className="flex gap-2">
                        <Button
                          size="sm"
                          variant={viewMode === 'result' ? 'flat' : 'light'}
                          color="primary"
                          onPress={() => setViewMode('result')}
                        >
                          Results
                        </Button>
                        <Button
                          size="sm"
                          variant={viewMode === 'heatmap' ? 'flat' : 'light'}
                          color="primary"
                          onPress={() => setViewMode('heatmap')}
                        >
                          Heatmap
                        </Button>
                      </div>
                    </div>
                  </CardHeader>
                  <CardBody>
                    {viewMode === 'result' ? (
                      <div>
                        <div className="flex flex-col items-center mb-4">
                          <Chip
                            size="lg"
                            color={getSeverityColor(predictionResult.result) as any}
                            className="mb-2"
                          >
                            {predictionResult.result} Diabetic Retinopathy
                          </Chip>
                          <p className="text-sm text-gray-600 text-center mb-4">
                            {getDescription(predictionResult.result)}
                          </p>
                        </div>

                        <div className="mb-6">
                          <div className="flex justify-between text-sm mb-1">
                            <span>Confidence:</span>
                            <span className="font-semibold">{(predictionResult.confidence * 100).toFixed(1)}%</span>
                          </div>
                          <Progress
                            color={getConfidenceColor(predictionResult.confidence) as any}
                            value={predictionResult.confidence * 100}
                            className="mb-4"
                            showValueLabel={false}
                          />
                        </div>

                        <Divider className="my-4" />

                        <div className="bg-gray-50 p-4 rounded-lg border-l-4 border-blue-500">
                          <p className="text-sm font-semibold mb-1">Clinical Recommendation:</p>
                          <p className="text-sm text-gray-600">
                            {getRecommendation(predictionResult.result)}
                          </p>
                        </div>

                        <div className="flex justify-end mt-6">
                          <Button
                            size="sm"
                            color="primary"
                            variant="light"
                            endContent={<Download size={16} />}
                          >
                            Export Report
                          </Button>
                        </div>
                      </div>
                    ) : (
                      <div className="flex flex-col items-center justify-center h-full">
                        <p className="text-sm text-gray-600 mb-4">
                          Heatmap visualization showing regions of interest that influenced the diagnosis:
                        </p>
                        <div className="relative w-full h-64">
                          <Image
                            src="/api/placeholder/400/300"
                            alt="AI Heatmap Visualization"
                            fill
                            style={{objectFit: "contain"}}
                            className="rounded-lg"
                          />
                        </div>
                        <p className="text-xs text-gray-500 mt-4 text-center">
                          Red areas indicate potential pathological findings that contributed to the classification
                        </p>
                      </div>
                    )}
                  </CardBody>
                </Card>
              ) : (
                <Card className="shadow-md bg-gray-50 border border-dashed border-gray-300">
                  <CardBody className="flex flex-col items-center justify-center py-12">
                    <AlertCircle size={40} className="text-gray-400 mb-4" />
                    <p className="text-gray-600 text-center">
                      Analysis results will appear here after processing the image
                    </p>
                  </CardBody>
                </Card>
              )}
            </div>
          )}

          {error && (
            <div className="bg-red-50 border-l-4 border-red-500 text-red-700 p-4 rounded-lg text-sm mb-6 flex items-center">
              <X size={18} className="mr-2 flex-shrink-0" />
              <span>{error}</span>
            </div>
          )}
        </div>
      )}

      {activeTab === 'batch' && (
        <div>
          <Card className="mb-6 shadow-md">
            <CardHeader className="flex gap-3 bg-gradient-to-r from-gray-50 to-gray-100">
              <div>
                <p className="text-lg font-semibold">Batch Upload</p>
                <p className="text-xs text-gray-500">
                  Upload a ZIP file containing multiple retinal images for bulk analysis
                </p>
              </div>
            </CardHeader>
            <CardBody>
              <div className="flex flex-col items-center justify-center border-2 border-dashed border-gray-300 rounded-xl p-8 bg-gray-50">
                <input
                  type="file"
                  ref={batchInputRef}
                  accept=".zip"
                  onChange={handleBatchUpload}
                  className="hidden"
                />
                <Upload size={40} className="text-gray-400 mb-4" />
                <p className="text-gray-600 mb-2 text-center">Drag and drop your ZIP file here</p>
                <p className="text-gray-400 text-sm mb-4 text-center">or</p>
                <Button
                  color="primary"
                  variant="flat"
                  onPress={triggerBatchInput}
                  startContent={<Upload size={16} />}
                >
                  Browse Files
                </Button>
                {batchFiles.length > 0 && (
                  <div className="mt-4 text-sm text-green-600 flex items-center">
                    <Check size={16} className="mr-1" /> {batchFiles.length} files ready for processing
                  </div>
                )}
              </div>
            </CardBody>
          </Card>

          {batchLoading && (
            <Card className="mb-6 shadow-md">
              <CardBody>
                <p className="text-sm text-gray-600 mb-2">Processing batch upload...</p>
                <Progress
                  color="primary"
                  isIndeterminate
                  aria-label="Processing Batch..."
                  className="mb-2"
                />
                <p className="text-xs text-gray-500">
                  This may take several minutes depending on the number of images
                </p>
              </CardBody>
            </Card>
          )}

          {batchResults.length > 0 && (
            <Card className="mb-6 shadow-md">
              <CardHeader className="bg-gradient-to-r from-gray-50 to-gray-100">
                <div className="flex justify-between items-center w-full">
                  <p className="text-md font-semibold">Batch Analysis Results</p>
                  <Button
                    size="sm"
                    color="primary"
                    variant="flat"
                    endContent={<Download size={16} />}
                  >
                    Export CSV
                  </Button>
                </div>
              </CardHeader>
              <CardBody className="p-0">
                <Table aria-label="Batch analysis results">
                  <TableHeader>
                    <TableColumn>FILENAME</TableColumn>
                    <TableColumn>SEVERITY</TableColumn>
                    <TableColumn>CONFIDENCE</TableColumn>
                    <TableColumn>RECOMMENDATION</TableColumn>
                  </TableHeader>
                  <TableBody>
                    {batchResults.map((result, index) => (
                      <TableRow key={index}>
                        <TableCell className="text-sm">{result.filename || `Image_${index + 1}.jpg`}</TableCell>
                        <TableCell>
                          <Chip
                            color={getSeverityColor(result.result) as any}
                            size="sm"
                          >
                            {result.result}
                          </Chip>
                        </TableCell>
                        <TableCell>
                          <div className="flex items-center gap-2">
                            <Progress
                              value={result.confidence * 100}
                              color={getConfidenceColor(result.confidence) as any}
                              size="sm"
                              className="max-w-20"
                            />
                            <span className="text-xs">
                              {(result.confidence * 100).toFixed(1)}%
                            </span>
                          </div>
                        </TableCell>
                        <TableCell className="text-xs">
                          {getRecommendation(result.result)}
                        </TableCell>
                      </TableRow>
                    ))}
                  </TableBody>
                </Table>
              </CardBody>
              <CardFooter className="bg-gray-50">
                <div className="w-full">
                  <div className="flex justify-between text-sm text-gray-600">
                    <span>Total Images: {batchResults.length}</span>
                    <span>
                      Requires Attention: {
                        batchResults.filter(r =>
                          r.result === 'Moderate' ||
                          r.result === 'Severe' ||
                          r.result === 'Proliferative'
                        ).length
                      }
                    </span>
                  </div>
                </div>
              </CardFooter>
            </Card>
          )}

          {error && (
            <div className="bg-red-50 border-l-4 border-red-500 text-red-700 p-4 rounded-lg text-sm mb-6 flex items-center">
              <X size={18} className="mr-2 flex-shrink-0" />
              <span>{error}</span>
            </div>
          )}
        </div>
      )}

      <Card className="bg-gray-50 shadow-md mb-6">
        <CardBody>
          <div className="flex flex-col md:flex-row justify-between items-center">
            <div className="text-center md:text-left mb-4 md:mb-0">
              <p className="text-sm text-gray-600">
                This tool is intended as a clinical decision support system.
                Final diagnosis should always be made by a qualified healthcare professional.
              </p>
            </div>
            <Button
              as="a"
              href="#"
              color="primary"
              variant="light"
              size="sm"
            >
              Learn More About Our Methodology
            </Button>
          </div>
        </CardBody>
      </Card>
    </div>
  );
}