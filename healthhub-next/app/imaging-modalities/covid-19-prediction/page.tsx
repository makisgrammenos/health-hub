'use client';

import React, { useState } from 'react';
import axios from 'axios';
import { Card, CardBody, CardFooter } from '@nextui-org/card';
import { Table, TableHeader, TableBody, TableColumn, TableRow, TableCell } from '@nextui-org/table';
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, PieChart,RadarChart, Pie, Cell, Legend } from 'recharts';
import Image from 'next/image';
import { Button } from '@nextui-org/button';

export default function Covid19Classification() {
  const [uploadedImage, setUploadedImage] = useState<string | null>(null);
  const [file, setFile] = useState<File | null>(null);
  const [prediction, setPrediction] = useState<string | null>(null);
  const [probabilities, setProbabilities] = useState<number[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleImageUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const selectedFile = event.target.files?.[0];
    if (!selectedFile) return;

    setError(null);
    setUploadedImage(URL.createObjectURL(selectedFile));
    setFile(selectedFile);
    setPrediction(null);
    setProbabilities([]);
  };

  const predictCovid19 = async () => {
    if (!file) {
      setError('Please upload an image.');
      return;
    }

    const formData = new FormData();
    formData.append('file', file);

    try {
      setLoading(true);
      const response = await axios.post('http://localhost:8000/imaging/covid/predict-covid', formData, {
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

  const transformProbabilitiesForChart = () => {
    const labels = ['Normal', 'Pneumonia', 'COVID-19'];
    return probabilities.map((prob, idx) => ({
      name: labels[idx],
      probability: prob * 100,
    }));
  };

  const COLORS = ['#0088FE', '#00C49F', '#FFBB28'];

  const topPrediction = prediction;

  return (
    <div className="max-w-7xl mx-auto p-8">
      <h1 className="text-4xl font-extrabold text-center mb-8">
        🦠 COVID-19 Classification
      </h1>
      <p className="text-lg text-gray-600 text-center mb-12">
        Use our AI-powered model to classify chest X-ray images and identify potential signs of COVID-19, Pneumonia, or Normal conditions.
      </p>

      <section>
        <h2 className="text-2xl font-bold mb-6">Upload a Chest X-Ray</h2>
        <input
          type="file"
          accept="image/*"
          onChange={handleImageUpload}
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
                onPress={predictCovid19}
                color="primary"
              >
                Predict
              </Button>
            </CardFooter>
          </Card>
        )}

        {prediction && probabilities.length > 0 && (
          <div>
            <h3 className="text-lg font-bold mb-4">Visualization</h3>
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-8 mb-6">
              {/* Bar Chart */}
              <div className="h-[500px] w-full">
                <ResponsiveContainer width="100%" height="100%">
                  <BarChart data={transformProbabilitiesForChart()} barCategoryGap="10%">
                    <CartesianGrid strokeDasharray="3 3" />
                    <XAxis dataKey="name" />
                    <YAxis />
                    <Tooltip />
                    <Bar
                      dataKey="probability"
                      fill="#8884d8"
                      label={{ position: 'top', fontSize: 12 }}
                    >
                      {transformProbabilitiesForChart().map((entry, index) => (
                        <Cell
                          key={`cell-${index}`}
                          fill={topPrediction === entry.name ? '#FF0000' : COLORS[index % COLORS.length]}
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
                      data={transformProbabilitiesForChart()}
                      dataKey="probability"
                      nameKey="name"
                      cx="50%"
                      cy="50%"
                      outerRadius={180}
                      innerRadius={100}
                      fill="#8884d8"
                      label={({ name, probability }) => `${name}: ${probability.toFixed(1)}%`}
                      labelLine={false}
                    >
                      {transformProbabilitiesForChart().map((entry, index) => (
                        <Cell
                          key={`cell-${index}`}
                          fill={topPrediction === entry.name ? '#FF0000' : COLORS[index % COLORS.length]}
                        />
                      ))}
                    </Pie>
                    <Tooltip />
                    <Legend verticalAlign="bottom" height={36} />
                  </PieChart>
                </ResponsiveContainer>
              </div>
              
            </div>

            {/* Prediction Results Table */}
            <Table aria-label="Prediction Results" className="mt-6">
              <TableHeader>
                <TableColumn>Condition</TableColumn>
                <TableColumn>Probability</TableColumn>
              </TableHeader>
              <TableBody>
                {transformProbabilitiesForChart().map((result, idx) => (
                  <TableRow
                    key={idx}
                    className={topPrediction === result.name ? 'bg-yellow-100 font-bold' : ''}
                  >
                    <TableCell>{result.name}</TableCell>
                    <TableCell>{result.probability.toFixed(2)}%</TableCell>
                  </TableRow>
                ))}
              </TableBody>
            </Table>
          </div>
        )}

        {error && <p className="text-red-600 mt-4 text-center">{error}</p>}
      </section>
    </div>
  );
}
