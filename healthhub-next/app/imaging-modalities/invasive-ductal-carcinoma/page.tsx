'use client';

import React, { useState } from 'react';
import axios from 'axios';
import { Tab, Tabs } from '@nextui-org/tabs';
import { Card, CardBody } from '@nextui-org/card';
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, } from 'recharts';
import { Button } from '@nextui-org/button';
import { Dropdown,DropdownItem } from '@nextui-org/react';

export default function BreastCancerPrediction() {
  const [uploadedImage, setUploadedImage] = useState<string | null>(null);
  const [file, setFile] = useState<File | null>(null);
  const [batchFile, setBatchFile] = useState<File | null>(null);
  const [prediction, setPrediction] = useState<string | null>(null);
  const [probabilities, setProbabilities] = useState<{ [key: string]: number }>({});
  const [batchResults, setBatchResults] = useState([]);
  const [selectedPatch, setSelectedPatch] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleImageUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const selectedFile = event.target.files?.[0];
    if (!selectedFile) return;

    setError(null);
    setFile(selectedFile);
    setUploadedImage(URL.createObjectURL(selectedFile));
    setPrediction(null);
    setProbabilities({});
  };

  const handleBatchUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const selectedFile = event.target.files?.[0];
    if (!selectedFile) return;

    setError(null);
    setBatchFile(selectedFile);
    setBatchResults([]);
    setSelectedPatch(null);
  };

  const predictBreastCancer = async () => {
    if (!file) {
      setError('Please upload an image.');
      return;
    }

    const formData = new FormData();
    formData.append('file', file);

    try {
      setLoading(true);
      const response = await axios.post('http://localhost:8000/imaging/breast-cancer/predict', formData, {
        headers: { 'Content-Type': 'multipart/form-data' },
      });

      setPrediction(response.data.result);
      setProbabilities(response.data.probabilities);
    } catch (err) {
      setError('Error occurred while predicting the image.');
    } finally {
      setLoading(false);
    }
  };

  const predictBatchBreastCancer = async () => {
    if (!batchFile) {
      setError('Please upload a ZIP file.');
      return;
    }

    const formData = new FormData();
    formData.append('file', batchFile);

    try {
      setLoading(true);
      const response = await axios.post('http://localhost:8000/imaging/breast-cancer/predict-zip', formData, {
        headers: { 'Content-Type': 'multipart/form-data' },
      });

      setBatchResults(response.data.results);
    } catch (err) {
      setError('Error occurred while processing the batch upload.');
    } finally {
      setLoading(false);
    }
  };

  const transformBatchResultsForChart = () => {
    const aggregatedResults = batchResults.reduce(
      (acc, result) => {
        if (result.probabilities && result.probabilities.length > 0) {
          const [noIDC, IDC] = result.probabilities[0];
          acc['No IDC'] += noIDC;
          acc['IDC'] += IDC;
        }
        return acc;
      },
      { 'No IDC': 0, 'IDC': 0 }
    );

    const totalPatches = batchResults.filter(result => result.probabilities).length;
    return Object.entries(aggregatedResults).map(([key, value]) => ({
      name: key,
      probability: (value / totalPatches) * 100,
    }));
  };

  const renderSelectedPatchChart = () => {
    if (!selectedPatch || !selectedPatch.probabilities) return null;

    const patchData = selectedPatch.probabilities[0].map((prob, idx) => ({
      name: idx === 0 ? 'No IDC' : 'IDC',
      probability: prob * 100,
    }));

    return (
      <ResponsiveContainer width="100%" height={300}>
        <BarChart data={patchData}>
          <CartesianGrid strokeDasharray="3 3" />
          <XAxis dataKey="name" />
          <YAxis />
          <Tooltip />
          <Bar dataKey="probability" fill="#8884d8" />
        </BarChart>
      </ResponsiveContainer>
    );
  };

  return (
    <div className="max-w-7xl mx-auto p-8">
      <h1 className="text-4xl font-extrabold text-center mb-8">🎗️ Breast Cancer Prediction</h1>

      <Tabs aria-label="Prediction Type" className="mb-6">
        <Tab title="Single Prediction">
          <section>
            <h2 className="text-2xl font-bold mb-6">Upload Breast Tissue Image</h2>
            <input
              type="file"
              accept="image/*"
              onChange={handleImageUpload}
              className="block w-full max-w-md mx-auto text-sm text-gray-500 file:mr-4 file:py-2 file:px-4 file:rounded-full file:border-0 file:bg-pink-50 file:text-pink-700 hover:file:bg-pink-100 mb-4"
            />
            <Button isDisabled={loading} isLoading={loading} onPress={predictBreastCancer} color="primary">
              Predict
            </Button>
            {uploadedImage && prediction && (
              <Card className="mt-6">
                <CardBody>
                  <p>Prediction: {prediction}</p>
                  {Object.entries(probabilities).map(([label, prob]) => (
                    <p key={label}>
                      {label}: {(prob * 100).toFixed(2)}%
                    </p>
                  ))}
                </CardBody>
              </Card>
            )}
          </section>
        </Tab>

        <Tab title="Batch Prediction">
          <section>
            <h2 className="text-2xl font-bold mb-6">Upload ZIP of Image Patches</h2>
            <input
              type="file"
              accept=".zip"
              onChange={handleBatchUpload}
              className="block w-full max-w-md mx-auto text-sm text-gray-500 file:mr-4 file:py-2 file:px-4 file:rounded-full file:border-0 file:bg-blue-50 file:text-blue-700 hover:file:bg-blue-100 mb-4"
            />
            <Button isDisabled={loading} isLoading={loading} onPress={predictBatchBreastCancer} color="primary">
              Predict Batch
            </Button>

            {batchResults.length > 0 && (
              <div className="mt-8">
                <h3 className="text-lg font-bold mb-4">Aggregated Visualization</h3>
                <ResponsiveContainer width="100%" height={400}>
                  <BarChart data={transformBatchResultsForChart()} barCategoryGap="10%">
                    <CartesianGrid strokeDasharray="3 3" />
                    <XAxis dataKey="name" />
                    <YAxis />
                    <Tooltip />
                    <Bar dataKey="probability" fill="#8884d8" />
                  </BarChart>
                </ResponsiveContainer>

                <h3 className="text-lg font-bold mb-4 mt-6">Select a Patch</h3>
                <Dropdown>
                  {batchResults.map((result, idx) => (
                    
                    <DropdownItem key={idx} onClick={() => setSelectedPatch(result)}>
                      {result.filename}
                    </DropdownItem>
                  ))}
                </Dropdown>

                {selectedPatch && (
                  <div className="mt-6">
                    <h3 className="text-lg font-bold">Patch Details: {selectedPatch.filename}</h3>
                    {renderSelectedPatchChart()}
                  </div>
                )}
              </div>
            )}
          </section>
        </Tab>
      </Tabs>

      {error && <p className="text-red-600 mt-4 text-center">{error}</p>}
    </div>
  );
}
