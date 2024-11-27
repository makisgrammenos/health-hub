'use client';

import React, { useState } from 'react';
import axios from 'axios';
import Image from 'next/image';
import { Card, CardBody, CardFooter } from '@nextui-org/card';
import { Button } from '@nextui-org/button';
import { Progress } from '@nextui-org/progress';
import { Tooltip } from '@nextui-org/tooltip';
import { Switch,cn } from '@nextui-org/react';

// Define types for responses and states
interface PredictionResult {
  result: string;
  confidence: number;
}

interface BatchPrediction {
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

  // Handle single image upload
  const handleImageUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const selectedFile = event.target.files?.[0];
    if (!selectedFile) return;

    setError(null);
    setUploadedImage(URL.createObjectURL(selectedFile));
    setFile(selectedFile);
    setPredictionResult(null); // Reset previous results
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
      setPredictionResult(response.data);
      console.log(response.data);
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
    } catch (err) {
      setError('An error occurred while processing the batch upload. Please try again.');
    } finally {
      setBatchLoading(false);
    }
  };

  return (
    <div className="max-w-7xl mx-auto p-8">
      <h1 className="text-3xl font-bold text-center mb-4">
        Diabetic Retinopathy Classification
      </h1>
      <p className="text-lg text-gray-600 text-center mb-6">
        Upload your retinal images to classify them for diabetic retinopathy using advanced deep
        learning models.
      </p>

      {/* Custom Switch to toggle between Single Image and Batch Upload */}
      <div className="flex justify-center mb-6">
        <Switch
          checked={activeTab === 'batch'}
          onChange={(e) => setActiveTab(e.target.checked ? 'batch' : 'single')}
          color="primary"
          size="lg"
          classNames={{
            base: cn(
              'inline-flex flex-row-reverse w-full max-w-md bg-content1 hover:bg-content2 items-center',
              'justify-between cursor-pointer rounded-lg gap-2 p-4 border-2 border-transparent',
              'data-[selected=true]:border-primary'
            ),
            wrapper: 'p-0 h-4 overflow-visible',
            thumb: cn(
              'w-6 h-6 border-2 shadow-lg',
              'group-data-[hover=true]:border-primary',
              // selected
              'group-data-[selected=true]:ml-6',
              // pressed
              'group-data-[pressed=true]:w-7',
              'group-data-[selected]:group-data-[pressed]:ml-4'
            ),
          }}
        >
          <div className="flex flex-col gap-1">
            <p className="text-medium">
              {activeTab === 'batch' ? 'Batch Upload' : 'Single Image'}
            </p>
            <p className="text-tiny text-default-400">
              {activeTab === 'batch'
                ? 'Upload multiple images in a ZIP file for batch classification.'
                : 'Upload a single retinal image for classification.'}
            </p>
          </div>
        </Switch>
      </div>

      {/* Section Heading */}
      <h2 className="text-2xl font-bold mb-6 text-center">
        {activeTab === 'single' ? 'Upload a Single Image' : 'Batch Upload (ZIP File)'}
      </h2>

      {activeTab === 'single' && (
        <div>
          <Card className="mb-6">
            <CardBody>
              <input
                type="file"
                accept="image/*"
                onChange={handleImageUpload}
                className="block w-full max-w-sm mx-auto text-sm text-gray-500
                file:mr-4 file:py-2 file:px-4 file:rounded-full file:border-0
                file:bg-indigo-600 file:text-white hover:file:bg-indigo-700 mb-4 cursor-pointer"
              />
            </CardBody>
          </Card>

          {uploadedImage && (
            <Card className="mb-6">
              <CardBody>
                <Image
                  src={uploadedImage}
                  alt="Uploaded Image"
                  width={600}
                  height={400}
                  className="object-cover rounded-lg mx-auto"
                />
              </CardBody>
              <CardFooter className="flex justify-center">
                <Button
                  isDisabled={loading}
                  isLoading={loading}
                  onPress={handlePrediction}
                  variant="shadow"
                  color="primary"
                >
                  {loading ? 'Processing...' : 'Predict'}
                </Button>
              </CardFooter>
            </Card>
          )}

          {predictionResult && (
            <Card className="mb-6">
              <CardBody>
                <h2 className="text-xl font-bold text-center mb-2">Prediction Results</h2>
                <p
                  className={`text-lg font-semibold text-center ${
                    predictionResult.result === 'Severe' ? 'text-red-600' : 'text-green-600'
                  }`}
                >
                  Prediction: {predictionResult.result}
                </p>
                <Progress
                  color="primary"
                  value={predictionResult.confidence * 100}
                  className="mt-4"
                />
              </CardBody>
            </Card>
          )}

          {error && (
            <div className="bg-red-100 text-red-700 px-4 py-3 rounded-lg text-center mb-6">
              {error}
            </div>
          )}
        </div>
      )}

      {activeTab === 'batch' && (
        <div>
          <Card className="mb-6">
            <CardBody>
              <input
                type="file"
                accept=".zip"
                onChange={handleBatchUpload}
                className="block w-full max-w-sm mx-auto text-sm text-gray-500
                file:mr-4 file:py-2 file:px-4 file:rounded-full file:border-0
                file:bg-green-600 file:text-white hover:file:bg-green-700 mb-4 cursor-pointer"
              />
            </CardBody>
          </Card>

          {batchLoading && (
            <Progress
              color="primary"
              value={100}
              label="Processing Batch..."
              className="mb-6"
            />
          )}

          {batchResults.length > 0 && (
            <Card className="mb-6">
              <CardBody>
                <h2 className="text-xl font-bold text-center mb-4">Batch Prediction Results</h2>
                <ul>
                  {batchResults.map((result, index) => (
                    <li key={index} className="text-lg mb-2">
                      <b>Image {index + 1}:</b> {result.result} (
                      {(result.confidence * 100).toFixed(2)}%)
                    </li>
                  ))}
                </ul>
              </CardBody>
            </Card>
          )}

          {error && (
            <div className="bg-red-100 text-red-700 px-4 py-3 rounded-lg text-center mb-6">
              {error}
            </div>
          )}
        </div>
      )}

      <Tooltip content="Sample Images Coming Soon" placement="bottom">
        <Button isDisabled className="mt-6 mx-auto block">
          View Sample Images
        </Button>
      </Tooltip>
    </div>
  );
}
