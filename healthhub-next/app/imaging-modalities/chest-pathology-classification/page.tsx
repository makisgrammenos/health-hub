'use client';

import React, { useState } from 'react';
import axios from 'axios';
import { Card, CardBody, CardFooter } from '@nextui-org/card';
import { Button } from '@nextui-org/button';
import { Progress } from '@nextui-org/progress';
import { Table } from '@nextui-org/table';
import Image from 'next/image';

export default function ChestXRayClassification() {
  const [activeTab, setActiveTab] = useState('single');
  const [uploadedImage, setUploadedImage] = useState<string | null>(null);
  const [file, setFile] = useState<File | null>(null);
  const [predictionResults, setPredictionResults] = useState([]);
  const [batchFiles, setBatchFiles] = useState<FileList | null>(null);
  const [batchResults, setBatchResults] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleSingleImageUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const selectedFile = event.target.files?.[0];
    if (!selectedFile) return;

    setError(null);
    setUploadedImage(URL.createObjectURL(selectedFile));
    setFile(selectedFile);
    setPredictionResults([]);
  };

  const handleBatchUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const files = event.target.files;
    if (!files || files.length === 0) return;

    setError(null);
    setBatchFiles(files);
    setBatchResults([]);
  };

  const predictSingleImage = async () => {
    if (!file) {
      setError('Please upload an image.');
      return;
    }

    const formData = new FormData();
    formData.append('file', file);

    try {
      setLoading(true);
      const response = await axios.post('/api/chest-xray/predict', formData, {
        headers: { 'Content-Type': 'multipart/form-data' },
      });
      setPredictionResults(response.data.results);
    } catch (err) {
      setError('Error occurred while predicting the image.');
    } finally {
      setLoading(false);
    }
  };

  const predictBatchImages = async () => {
    if (!batchFiles) {
      setError('Please upload a ZIP file.');
      return;
    }

    const formData = new FormData();
    Array.from(batchFiles).forEach((file) => formData.append('files', file));

    try {
      setLoading(true);
      const response = await axios.post('/api/chest-xray/batch-predict', formData, {
        headers: { 'Content-Type': 'multipart/form-data' },
      });
      setBatchResults(response.data.predictions);
    } catch (err) {
      setError('Error occurred while processing the batch upload.');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="max-w-7xl mx-auto p-8">
      <h1 className="text-4xl font-extrabold text-center mb-8">
        🩺 Chest X-Ray Classification
      </h1>
      <p className="text-lg text-gray-600 text-center mb-12">
        Use our AI-powered model to classify chest X-ray images and identify potential pathologies.
      </p>

      <div className="flex justify-center space-x-4 mb-8">
        <Button
          color={activeTab === 'single' ? 'primary' : 'default'}
          onPress={() => setActiveTab('single')}
        >
          Single Image
        </Button>
        <Button
          color={activeTab === 'batch' ? 'primary' : 'default'}
          onPress={() => setActiveTab('batch')}
        >
          Batch Prediction
        </Button>
      </div>

      {activeTab === 'single' && (
        <section>
          <h2 className="text-2xl font-bold mb-6">Upload a Chest X-Ray</h2>
          <input
            type="file"
            accept="image/*, .dcm"
            onChange={handleSingleImageUpload}
            className="block w-full max-w-md mx-auto text-sm text-gray-500 file:mr-4 file:py-2 file:px-4 file:rounded-full file:border-0 file:bg-indigo-50 file:text-indigo-700 hover:file:bg-indigo-100 mb-4"
          />
          {uploadedImage && (
            <Card className="mb-6">
              <CardBody>
                <Image
                  src={uploadedImage}
                  alt="Uploaded X-Ray"
                  width={600}
                  height={400}
                  className="rounded-lg"
                />
              </CardBody>
              <CardFooter className="flex justify-center">
                <Button
                  isDisabled={loading}
                  isLoading={loading}
                  onPress={predictSingleImage}
                  color="primary"
                >
                  Predict
                </Button>
              </CardFooter>
            </Card>
          )}

          {predictionResults.length > 0 && (
            <Table
              aria-label="Prediction Results"
              className="mt-6"
            >
              <thead>
                <tr>
                  <th>Pathology</th>
                  <th>Probability</th>
                </tr>
              </thead>
              <tbody>
                {predictionResults.map((result: any, idx: number) => (
                  <tr key={idx}>
                    <td>{result.Pathology}</td>
                    <td>{result.Probability}</td>
                  </tr>
                ))}
              </tbody>
            </Table>
          )}

          {error && <p className="text-red-600 mt-4 text-center">{error}</p>}
        </section>
      )}

      {activeTab === 'batch' && (
        <section>
          <h2 className="text-2xl font-bold mb-6">Batch Prediction</h2>
          <input
            type="file"
            accept=".zip"
            onChange={handleBatchUpload}
            className="block w-full max-w-md mx-auto text-sm text-gray-500 file:mr-4 file:py-2 file:px-4 file:rounded-full file:border-0 file:bg-indigo-50 file:text-indigo-700 hover:file:bg-indigo-100 mb-4"
          />
          {batchFiles && (
            <Button
              isDisabled={loading}
              isLoading={loading}
              onPress={predictBatchImages}
              color="primary"
              className="block mx-auto"
            >
              Predict Batch
            </Button>
          )}

          {batchResults.length > 0 && (
            <Table
              aria-label="Batch Predictions"
              className="mt-6"
            >
              <thead>
                <tr>
                  <th>Image</th>
                  <th>Top Pathology</th>
                  <th>Probability</th>
                </tr>
              </thead>
              <tbody>
                {batchResults.map((result: any, idx: number) => (
                  <tr key={idx}>
                    <td>{result.Image}</td>
                    <td>{result['Top Pathology']}</td>
                    <td>{result.Probability}</td>
                  </tr>
                ))}
              </tbody>
            </Table>
          )}

          {error && <p className="text-red-600 mt-4 text-center">{error}</p>}
        </section>
      )}
    </div>
  );
}
