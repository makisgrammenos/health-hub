'use client';

import React, { useState } from 'react';
import axios from 'axios';
import { Dropdown } from '@nextui-org/react';
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts';
import { Button } from '@nextui-org/button';

export default function SkinCancerPrediction() {
  const [uploadedFile, setUploadedFile] = useState<File | null>(null);
  const [selectedFile, setSelectedFile] = useState<string | null>(null);
  const [prediction, setPrediction] = useState<string | null>(null);
  const [probabilities, setProbabilities] = useState<{ [key: string]: number } | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;

    setError(null);
    setUploadedFile(file);
    setSelectedFile(file.name);
    setPrediction(null);
    setProbabilities(null);
  };

  const handlePrediction = async () => {
    if (!uploadedFile) {
      setError('Please upload an image.');
      return;
    }

    const formData = new FormData();
    formData.append('file', uploadedFile);

    try {
      setLoading(true);
      const response = await axios.post('http://localhost:8000/imaging/skin-cancer/predict', formData, {
        headers: { 'Content-Type': 'multipart/form-data' },
      });

      setPrediction(response.data.prediction);
      setProbabilities(response.data.probabilities);
    } catch (err) {
      setError('Error occurred while predicting the image.');
    } finally {
      setLoading(false);
    }
  };

  const chartData = probabilities
    ? Object.entries(probabilities).map(([key, value]) => ({ name: key, probability: value * 100 }))
    : [];

  return (
    <div className="max-w-7xl mx-auto p-8">
      <h1 className="text-4xl font-extrabold text-center mb-8">
        🖼️ Skin Cancer Prediction
      </h1>

      <section>
        <h2 className="text-2xl font-bold mb-6">Upload Skin Lesion Image</h2>
        <input
          type="file"
          accept="image/*"
          onChange={handleFileUpload}
          className="block w-full max-w-md mx-auto text-sm text-gray-500 file:mr-4 file:py-2 file:px-4 file:rounded-full file:border-0 file:bg-blue-50 file:text-blue-700 hover:file:bg-blue-100 mb-4"
        />
        <Button isDisabled={loading} isLoading={loading} onPress={handlePrediction} color="primary">
          Predict
        </Button>

        {selectedFile && <p className="mt-4">Selected File: {selectedFile}</p>}

        {prediction && probabilities && (
          <div className="mt-8">
            <h3 className="text-lg font-bold mb-4">Prediction Result</h3>
            <p className="text-lg">Prediction: {prediction}</p>

            <ResponsiveContainer width="100%" height={400}>
              <BarChart data={chartData} barCategoryGap="10%">
                <CartesianGrid strokeDasharray="3 3" />
                <XAxis dataKey="name" />
                <YAxis />
                <Tooltip />
                <Bar dataKey="probability" fill="#8884d8" />
              </BarChart>
            </ResponsiveContainer>
          </div>
        )}

        {error && <p className="text-red-600 mt-4 text-center">{error}</p>}
      </section>
    </div>
  );
}
