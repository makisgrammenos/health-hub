'use client';

import React, { useState } from 'react';
import axios from 'axios';
import { Card, CardBody, CardFooter } from '@nextui-org/card';
import { Table, TableHeader, TableBody, TableColumn, TableRow, TableCell } from '@nextui-org/table';
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, PieChart, Pie, Cell, Legend } from 'recharts';
import Image from 'next/image';
import { Button } from '@nextui-org/button';
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
      const response = await axios.post('http://localhost:8000/imaging/chest-x-ray/predict', formData, {
        headers: { 'Content-Type': 'multipart/form-data' },
      });
      setPredictionResults(response.data.predictions);
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
      const response = await axios.post('http://localhost:8000/imaging/chest-x-ray/predict', formData, {
        headers: { 'Content-Type': 'multipart/form-data' },
      });
      setBatchResults(response.data.predictions);
    } catch (err) {
      setError('Error occurred while processing the batch upload.');
    } finally {
      setLoading(false);
    }
  };

  const transformPredictionDataForChart = () => {
    return predictionResults.map((result: any) => ({
      name: result.Pathology,
      probability: parseFloat(result.Probability.replace('%', '')),
    }));
  };

  const getTopPrediction = () => {
    if (predictionResults.length === 0) return null;
    return predictionResults.reduce((top, current) =>
      parseFloat(current.Probability.replace('%', '')) > parseFloat(top.Probability.replace('%', '')) ? current : top
    );
  };

  const COLORS = [
    '#0088FE', '#00C49F', '#FFBB28', '#FF8042', '#AF19FF', '#FF6361', '#FFD700', '#FF4500', '#32CD32', '#6A5ACD',
  ];

  const topPrediction = getTopPrediction();

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
            <div>
              <h3 className="text-lg font-bold mb-4">Visualization</h3>
              <div className="grid grid-cols-1 lg:grid-cols-2 gap-8 mb-6">
  {/* Bar Chart */}
  <div className="h-[500px] w-full">
    <ResponsiveContainer width="100%" height="100%">
      <BarChart data={transformPredictionDataForChart()} barCategoryGap="10%">
        <CartesianGrid strokeDasharray="3 3" />
        <XAxis dataKey="name" />
        <YAxis />
        <Tooltip />
        <Bar
          dataKey="probability"
          fill="#8884d8"
          label={{ position: 'top', fontSize: 12 }}
        >
          {transformPredictionDataForChart().map((entry, index) => (
            <Cell
              key={`cell-${index}`}
              fill={
                topPrediction && entry.name === topPrediction.Pathology
                  ? '#FF0000'
                  : COLORS[index % COLORS.length]
              }
            />
          ))}
        </Bar>
      </BarChart>
    </ResponsiveContainer>
  </div>

  {/* Donut Chart */}
  <div className="h-[500px] w-full">
    <ResponsiveContainer width="100%" height="100%">
      <PieChart>
        <Pie
          data={transformPredictionDataForChart()}
          dataKey="probability"
          nameKey="name"
          cx="50%"
          cy="50%"
          outerRadius={180}
          innerRadius={100} // Donut effect
          fill="#8884d8"
          label={({ name, probability }) => `${name}: ${probability.toFixed(1)}%`}
          labelLine={false}
        >
          {transformPredictionDataForChart().map((entry, index) => (
            <Cell
              key={`cell-${index}`}
              fill={
                topPrediction && entry.name === topPrediction.Pathology
                  ? '#FF0000'
                  : COLORS[index % COLORS.length]
              }
            />
          ))}
        </Pie>
        <Tooltip />
        <Legend verticalAlign="bottom" height={36} />
      </PieChart>
    </ResponsiveContainer>
  </div>
</div>


              <Table aria-label="Prediction Results" className="mt-6">
                <TableHeader>
                  <TableColumn>Pathology</TableColumn>
                  <TableColumn>Probability</TableColumn>
                </TableHeader>
                <TableBody>
                  {predictionResults.map((result: any, idx: number) => (
                    <TableRow
                      key={idx}
                      className={
                        topPrediction && result.Pathology === topPrediction.Pathology
                          ? 'bg-yellow-100 font-bold'
                          : ''
                      }
                    >
                      <TableCell>{result.Pathology}</TableCell>
                      <TableCell>{result.Probability}</TableCell>
                    </TableRow>
                  ))}
                </TableBody>
              </Table>
            </div>
          )}

          {error && <p className="text-red-600 mt-4 text-center">{error}</p>}
        </section>
      )}
    </div>
  );
}
