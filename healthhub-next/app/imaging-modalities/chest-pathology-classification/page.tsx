'use client';

import React, { useState } from 'react';
import axios from 'axios';
import { Card, CardBody, CardFooter, CardHeader } from '@nextui-org/card';
import {
  Table,
  TableHeader,
  TableBody,
  TableColumn,
  TableRow,
  TableCell,
} from '@nextui-org/table';
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  PieChart,
  Pie,
  Cell,
  Legend,
  LabelList,
} from 'recharts';
import Image from 'next/image';
import { Button } from '@nextui-org/button';
import { Switch } from '@nextui-org/switch';
import { cn } from "@nextui-org/react";
import { Tabs, Tab } from "@nextui-org/tabs";
import { Progress } from "@nextui-org/progress";
import { Divider } from "@nextui-org/divider";
import { Chip } from "@nextui-org/chip";

export default function ChestXRayClassification() {
  const [uploadedImage, setUploadedImage] = useState(null);
  const [file, setFile] = useState(null);
  const [isDicomFile, setIsDicomFile] = useState(false);
  const [predictionResults, setPredictionResults] = useState([]);
  const [batchFiles, setBatchFiles] = useState(null);
  const [batchResults, setBatchResults] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [activeTab, setActiveTab] = useState('single');

  const handleSingleImageUpload = (event) => {
    const selectedFile = event.target.files?.[0];
    if (!selectedFile) return;

    setError(null);
    setFile(selectedFile);
    setPredictionResults([]);
    const fileExtension = selectedFile.name.split('.').pop()?.toLowerCase();

    if (fileExtension === 'dcm') {
      setIsDicomFile(true);
      setUploadedImage(null);
    } else {
      setIsDicomFile(false);
      setUploadedImage(URL.createObjectURL(selectedFile));
    }
  };

  const handleBatchUpload = (event) => {
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
      const response = await axios.post(
        'http://localhost:8000/imaging/chest-x-ray/predict',
        formData,
        {
          headers: { 'Content-Type': 'multipart/form-data' },
        }
      );
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
      const response = await axios.post(
        'http://localhost:8000/imaging/chest-x-ray/predict',
        formData,
        {
          headers: { 'Content-Type': 'multipart/form-data' },
        }
      );
      setBatchResults(response.data.predictions);
    } catch (err) {
      setError('Error occurred while processing the batch upload.');
    } finally {
      setLoading(false);
    }
  };

  const transformPredictionDataForChart = () => {
    return predictionResults.map((result) => ({
      name: result.Pathology,
      probability: parseFloat(result.Probability.replace('%', '')),
    }));
  };

  const getTopPrediction = () => {
    if (predictionResults.length === 0) return null;
    return predictionResults.reduce((top, current) =>
      parseFloat(current.Probability.replace('%', '')) >
      parseFloat(top.Probability.replace('%', ''))
        ? current
        : top
    );
  };

  const COLORS = [
    '#0088FE',
    '#00C49F',
    '#FFBB28',
    '#FF8042',
    '#AF19FF',
    '#FF6361',
    '#FFD700',
    '#FF4500',
    '#32CD32',
    '#6A5ACD',
  ];

  const topPrediction = getTopPrediction();

  // Calculate confidence level class based on highest prediction probability
  const getConfidenceClass = () => {
    if (!topPrediction) return "";
    const probability = parseFloat(topPrediction.Probability.replace('%', ''));
    if (probability >= 85) return "text-green-600";
    if (probability >= 70) return "text-yellow-600";
    return "text-red-600";
  };

  const getConfidenceText = () => {
    if (!topPrediction) return "";
    const probability = parseFloat(topPrediction.Probability.replace('%', ''));
    if (probability >= 85) return "High confidence";
    if (probability >= 70) return "Moderate confidence";
    return "Low confidence";
  };

  const formatDate = () => {
    const date = new Date();
    return new Intl.DateTimeFormat('en-US', {
      year: 'numeric',
      month: 'long',
      day: 'numeric',
      hour: 'numeric',
      minute: 'numeric',
      hour12: true
    }).format(date);
  };

  return (
    <div className="min-h-screen bg-gray-50">
      <header className="bg-blue-900 text-white py-4 px-8 shadow-md">
        <div className="max-w-7xl mx-auto flex justify-between items-center">
          <div className="flex items-center space-x-2">
            <svg xmlns="http://www.w3.org/2000/svg" className="h-8 w-8" viewBox="0 0 20 20" fill="currentColor">
              <path fillRule="evenodd" d="M4 5a2 2 0 00-2 2v8a2 2 0 002 2h12a2 2 0 002-2V7a2 2 0 00-2-2h-1.586a1 1 0 01-.707-.293l-1.121-1.121A2 2 0 0011.172 3H8.828a2 2 0 00-1.414.586L6.293 4.707A1 1 0 015.586 5H4zm6 9a3 3 0 100-6 3 3 0 000 6z" clipRule="evenodd" />
            </svg>
            <h1 className="text-2xl font-bold tracking-tight">Medical Imaging Analysis Platform</h1>
          </div>
          <div className="hidden md:flex items-center space-x-4">
            <span className="text-sm">{formatDate()}</span>
            <Chip color="primary" variant="flat">v2.4.1</Chip>
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto p-6 md:p-8">
        <Card className="shadow-lg border-0 mb-8">
          <CardHeader className="bg-gradient-to-r from-blue-800 to-blue-600 text-white">
            <div className="flex flex-col">
              <h2 className="text-2xl font-bold">Chest X-Ray Classification</h2>
              <p className="text-blue-100 text-sm mt-1">AI-powered analysis for radiological pathology detection</p>
            </div>
          </CardHeader>
          <CardBody>
            <div className="mb-6">
              <Tabs
                variant="bordered"
                selectedKey={activeTab}
                onSelectionChange={setActiveTab}
                className="w-full"
                aria-label="Analysis Options"
              >
                <Tab
                  key="single"
                  title={
                    <div className="flex items-center gap-2">
                      <svg xmlns="http://www.w3.org/2000/svg" className="h-4 w-4" viewBox="0 0 20 20" fill="currentColor">
                        <path d="M5 4a1 1 0 00-2 0v7.268a2 2 0 000 3.464V16a1 1 0 102 0v-1.268a2 2 0 000-3.464V4zM11 4a1 1 0 10-2 0v1.268a2 2 0 000 3.464V16a1 1 0 102 0V8.732a2 2 0 000-3.464V4zM16 3a1 1 0 011 1v7.268a2 2 0 010 3.464V16a1 1 0 11-2 0v-1.268a2 2 0 010-3.464V4a1 1 0 011-1z" />
                      </svg>
                      <span>Single Image Analysis</span>
                    </div>
                  }
                >
                  <div className="py-4">
                    <div className="max-w-7xl mx-auto">
                      <div className="bg-blue-50 rounded-lg p-4 mb-6">
                        <div className="flex items-start">
                          <div className="bg-blue-100 rounded-full p-2 mr-3">
                            <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5 text-blue-700" viewBox="0 0 20 20" fill="currentColor">
                              <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2h-1V9a1 1 0 00-1-1z" clipRule="evenodd" />
                            </svg>
                          </div>
                          <div>
                            <h3 className="font-semibold text-blue-800">Instructions</h3>
                            <p className="text-sm text-gray-600 mt-1">
                              Upload a single chest X-ray image for pathology prediction. The AI will analyze the image and identify potential conditions.
                            </p>
                            <p className="text-xs text-gray-500 mt-2">
                              Supported formats: JPEG, PNG, or DICOM (.dcm)
                            </p>
                          </div>
                        </div>
                      </div>

                      <div className="flex justify-center mb-4">
                        <label className="flex flex-col items-center w-full max-w-md px-4 py-6 bg-white text-blue-900 rounded-lg border-2 border-dashed border-blue-300 cursor-pointer hover:bg-blue-50 transition-colors">
                          <svg xmlns="http://www.w3.org/2000/svg" className="h-12 w-12 text-blue-400" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M7 16a4 4 0 01-.88-7.903A5 5 0 1115.9 6L16 6a5 5 0 011 9.9M15 13l-3-3m0 0l-3 3m3-3v12" />
                          </svg>
                          <span className="mt-2 text-base font-medium">Upload X-Ray Image</span>
                          <span className="mt-1 text-xs text-gray-500">
                            Click to browse or drag and drop
                          </span>
                          <input
                            type="file"
                            accept="image/*, .dcm"
                            onChange={handleSingleImageUpload}
                            className="hidden"
                          />
                        </label>
                      </div>

                      {(uploadedImage || file) && (
                        <Card className="mb-6 border shadow-sm">
                          <CardHeader className="bg-gray-50 border-b">
                            <div className="flex justify-between items-center w-full">
                              <h3 className="text-md font-semibold">Image Preview</h3>
                              <Chip size="sm" variant="flat">{file?.name}</Chip>
                            </div>
                          </CardHeader>
                          <CardBody>
                            {isDicomFile ? (
                              <div className="text-center py-12 bg-gray-50 rounded-lg">
                                <svg xmlns="http://www.w3.org/2000/svg" className="h-16 w-16 mx-auto text-gray-400" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
                                </svg>
                                <p className="text-md font-medium mt-4">DICOM file loaded successfully</p>
                                <p className="text-gray-500 mt-1 text-sm">
                                  File will be analyzed directly without preview
                                </p>
                              </div>
                            ) : (
                              <div className="flex justify-center border rounded-lg overflow-hidden">
                                <Image
                                  src={uploadedImage}
                                  alt="Uploaded X-Ray"
                                  width={600}
                                  height={400}
                                  className="w-full h-full object-contain max-h-[400px]"
                                />
                              </div>
                            )}
                          </CardBody>
                          <CardFooter className="flex justify-center border-t bg-gray-50">
                            <Button
                              isDisabled={loading}
                              isLoading={loading}
                              onPress={predictSingleImage}
                              color="primary"
                              className="px-8 font-medium"
                              startContent={
                                <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5" viewBox="0 0 20 20" fill="currentColor">
                                  <path fillRule="evenodd" d="M10 2a1 1 0 011 1v1.323l3.954 1.582 1.599-.8a1 1 0 01.894 1.79l-1.233.616 1.738 5.42a1 1 0 01-.285 1.05A3.989 3.989 0 0115 15a3.989 3.989 0 01-2.667-1.019 1 1 0 01-.285-1.05l1.715-5.349L11 6.477V16h2a1 1 0 110 2H7a1 1 0 110-2h2V6.477L6.237 7.582l1.715 5.349a1 1 0 01-.285 1.05A3.989 3.989 0 015 15a3.989 3.989 0 01-2.667-1.019 1 1 0 01-.285-1.05l1.738-5.42-1.233-.617a1 1 0 01.894-1.788l1.599.799L9 4.323V3a1 1 0 011-1z" clipRule="evenodd" />
                                </svg>
                              }
                            >
                              Analyze Image
                            </Button>
                          </CardFooter>
                        </Card>
                      )}

                    {predictionResults.length > 0 && (
                     <Card className="mt-8 border shadow-md w-full max-w-none max">

                          <CardHeader className="bg-gradient-to-r from-green-800 to-green-600 text-white">
                            <div className="flex flex-col w-full">
                              <div className="flex justify-between items-center">
                                <h3 className="text-xl font-bold">Analysis Results</h3>
<Chip
       variant="solid"
       color={
        // pick a solid color variant based on your same thresholds:
         topPrediction && parseFloat(topPrediction.Probability.replace('%','')) >= 85      ? 'success'
          : topPrediction && parseFloat(topPrediction.Probability.replace('%','')) >= 70
             ? 'warning'
            : 'danger'
       }
    >                                  {getConfidenceText()}
                                </Chip>
                              </div>
                              {topPrediction && (
                                <div className="mt-2 text-sm bg-black/20 rounded-md px-3 py-2">
                                  <span className="opacity-90">Primary finding:</span>{" "}
                                  <span className="font-bold">{topPrediction.Pathology}</span>
                                  <span className="ml-2">({topPrediction.Probability})</span>
                                </div>
                              )}
                            </div>
                          </CardHeader>
                          <CardBody>
                            <div className="grid grid-cols-1 lg:grid-cols-2 gap-8 mb-6">
                              {/* Bar Chart */}
                              <Card className="p-4 border shadow-sm">
                                <CardHeader className="px-0 pt-0 pb-2">
                                  <h4 className="text-md font-semibold">Probability Distribution (Bar Chart)</h4>
                                </CardHeader>
                                <CardBody className="px-0 py-0" >
                                  <div className="h-96">
                                    <ResponsiveContainer >
                                      <BarChart data={transformPredictionDataForChart()}>
                                        <CartesianGrid strokeDasharray="3 3" opacity={0.3} />
                                        <XAxis dataKey="name" tick={{ fontSize: 11 }} />
                                        <YAxis domain={[0, 100]} tickFormatter={(value) => `${value}%`} />
                                        <Tooltip formatter={(value) => [`${value.toFixed(1)}%`, 'Probability']} />
                                        <Bar dataKey="probability" radius={[4, 4, 0, 0]} barSize={60}>
                                          {transformPredictionDataForChart().map((entry, index) => (
                                            <Cell
                                              key={`cell-${index}`}
                                              fill={
                                                topPrediction &&
                                                entry.name === topPrediction.Pathology
                                                  ? '#16a34a'
                                                  : COLORS[index % COLORS.length]
                                              }
                                            />
                                          ))}
                                          <LabelList
                                            dataKey="probability"
                                            position="top"
                                            formatter={(value) => `${value.toFixed(1)}%`}
                                            style={{ fontSize: '11px', fill: '#374151' }}
                                          />
                                        </Bar>
                                      </BarChart>
                                    </ResponsiveContainer>
                                  </div>
                                </CardBody>
                              </Card>

                              {/* Donut Chart */}
                              <Card className="p-4 border shadow-sm">
                                <CardHeader className="px-0 pt-0 pb-2">
                                  <h4 className="text-md font-semibold">Findings Distribution (Pie Chart)</h4>
                                </CardHeader>
                                <CardBody className="px-0 py-0">
                                  <div className="h-72">
                                    <ResponsiveContainer>
                                      <PieChart>
                                        <Pie
                                          data={transformPredictionDataForChart()}
                                          dataKey="probability"
                                          nameKey="name"
                                          cx="50%"
                                          cy="50%"
                                          outerRadius={100}
                                          innerRadius={60}
                                          fill="#8884d8"
                                          labelLine={true}
                                          label={({ name, probability }) =>
                                            `${name}: ${probability.toFixed(1)}%`
                                          }
                                        >
                                          {transformPredictionDataForChart().map((entry, index) => (
                                            <Cell
                                              key={`cell-${index}`}
                                              fill={
                                                topPrediction &&
                                                entry.name === topPrediction.Pathology
                                                  ? '#16a34a'
                                                  : COLORS[index % COLORS.length]
                                              }
                                              stroke="#fff"
                                              strokeWidth={1}
                                            />
                                          ))}
                                        </Pie>
                                        <Tooltip formatter={(value) => [`${value.toFixed(1)}%`, 'Probability']} />
                                        <Legend verticalAlign="bottom" height={36} />
                                      </PieChart>
                                    </ResponsiveContainer>
                                  </div>
                                </CardBody>
                              </Card>
                            </div>

                            <Divider className="my-6" />

                            <div className="overflow-hidden">
                              <h4 className="text-md font-semibold mb-4">Detailed Analysis</h4>
                              <div className="bg-gray-50 p-4 rounded-lg mb-6">
                                {predictionResults.map((result, idx) => (
                                  <div key={idx} className="mb-3 last:mb-0">
                                    <div className="flex justify-between items-center mb-1">
                                      <div className="flex items-center">
                                        <span className={`font-medium ${topPrediction && result.Pathology === topPrediction.Pathology ? "text-green-700" : "text-gray-700"}`}>
                                          {result.Pathology}
                                        </span>
                                        {topPrediction && result.Pathology === topPrediction.Pathology && (
                                          <Chip size="sm" color="success" variant="flat" className="ml-2">Primary</Chip>
                                        )}
                                      </div>
                                      <span className="text-sm font-semibold">{result.Probability}</span>
                                    </div>
                                    <Progress
                                      value={parseFloat(result.Probability.replace('%', ''))}
                                      color={topPrediction && result.Pathology === topPrediction.Pathology ? "success" : "primary"}
                                      size="sm"
                                      radius="sm"
                                      classNames={{
                                        track: "bg-gray-200"
                                      }}
                                      aria-label={`${result.Pathology} probability`}
                                    />
                                  </div>
                                ))}
                              </div>

                              <Table
                                aria-label="Detailed Prediction Results"
                                className="shadow-sm"
                                selectionMode="none"
                                isHeaderSticky
                              >
                                <TableHeader>
                                  <TableColumn>Rank</TableColumn>
                                  <TableColumn>Pathology</TableColumn>
                                  <TableColumn>Classification</TableColumn>
                                  <TableColumn>Probability</TableColumn>
                                </TableHeader>
                                <TableBody>
                                  {predictionResults
                                    .sort((a, b) =>
                                      parseFloat(b.Probability.replace('%', '')) -
                                      parseFloat(a.Probability.replace('%', ''))
                                    )
                                    .map((result, idx) => {
                                      const probability = parseFloat(result.Probability.replace('%', ''));
                                      let severity = "Normal";
                                      let chipColor = "default";

                                      if (probability >= 80) {
                                        severity = "Critical";
                                        chipColor = "danger";
                                      } else if (probability >= 60) {
                                        severity = "Significant";
                                        chipColor = "warning";
                                      } else if (probability >= 40) {
                                        severity = "Moderate";
                                        chipColor = "primary";
                                      } else if (probability >= 20) {
                                        severity = "Mild";
                                        chipColor = "success";
                                      }

                                      return (
                                        <TableRow
                                          key={idx}
                                          className={
                                            topPrediction &&
                                            result.Pathology === topPrediction.Pathology
                                              ? 'bg-green-50'
                                              : ''
                                          }
                                        >
                                          <TableCell>{idx + 1}</TableCell>
                                          <TableCell>
                                            <div className="flex items-center">
                                              {result.Pathology}
                                              {idx === 0 && (
                                                <Chip size="sm" color="success" variant="dot" className="ml-2">Top</Chip>
                                              )}
                                            </div>
                                          </TableCell>
                                          <TableCell>
                                            <Chip size="sm" color={chipColor} variant="flat">{severity}</Chip>
                                          </TableCell>
                                          <TableCell>
                                            <div className="font-medium">
                                              {result.Probability}
                                            </div>
                                          </TableCell>
                                        </TableRow>
                                      );
                                    })}
                                </TableBody>
                              </Table>
                            </div>
                          </CardBody>
                          <CardFooter className="flex justify-center border-t pt-4 gap-2">
                            <Button
                              color="default"
                              variant="flat"
                              startContent={
                                <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5" viewBox="0 0 20 20" fill="currentColor">
                                  <path fillRule="evenodd" d="M6 2a2 2 0 00-2 2v12a2 2 0 002 2h8a2 2 0 002-2V7.414A2 2 0 0015.414 6L12 2.586A2 2 0 0010.586 2H6zm5 6a1 1 0 10-2 0v2H7a1 1 0 100 2h2v2a1 1 0 102 0v-2h2a1 1 0 100-2h-2V8z" clipRule="evenodd" />
                                </svg>
                              }
                            >
                              Export Report
                            </Button>
                            <Button
                              color="primary"
                              variant="flat"
                              startContent={
                                <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5" viewBox="0 0 20 20" fill="currentColor">
                                  <path d="M13.586 3.586a2 2 0 112.828 2.828l-.793.793-2.828-2.828.793-.793zM11.379 5.793L3 14.172V17h2.828l8.38-8.379-2.83-2.828z" />
                                </svg>
                              }
                            >
                              Add Notes
                            </Button>
                          </CardFooter>
                        </Card>
                      )}
                    </div>
                  </div>
                </Tab>
                <Tab
                  key="batch"
                  title={
                    <div className="flex items-center gap-2">
                      <svg xmlns="http://www.w3.org/2000/svg" className="h-4 w-4" viewBox="0 0 20 20" fill="currentColor">
                        <path d="M7 3a1 1 0 000 2h6a1 1 0 100-2H7zM4 7a1 1 0 011-1h10a1 1 0 110 2H5a1 1 0 01-1-1zM2 11a2 2 0 012-2h12a2 2 0 012 2v4a2 2 0 01-2 2H4a2 2 0 01-2-2v-4z" />
                      </svg>
                      <span>Batch Analysis</span>
                    </div>
                  }
                >
                  <div className="py-4">
                    <div className="max-w-xl mx-auto">
                      <div className="bg-blue-50 rounded-lg p-4 mb-6">
                        <div className="flex items-start">
                          <div className="bg-blue-100 rounded-full p-2 mr-3">
                            <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5 text-blue-700" viewBox="0 0 20 20" fill="currentColor">
                              <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2h-1V9a1 1 0 00-1-1z" clipRule="evenodd" />
                            </svg>
                          </div>
                          <div>
                            <h3 className="font-semibold text-blue-800">Batch Analysis Instructions</h3>
                            <p className="text-sm text-gray-600 mt-1">
                              Process multiple chest X-ray images at once. Upload a ZIP file containing X-ray images for efficient bulk analysis.
                            </p>
                            <p className="text-xs text-gray-500 mt-2">
                              Supported format: ZIP archive containing JPEG, PNG, or DICOM files
                            </p>
                          </div>
                        </div>
                      </div>

                      <div className="flex justify-center mb-4">
                        <label className="flex flex-col items-center w-full max-w-md px-4 py-6 bg-white text-blue-900 rounded-lg border-2 border-dashed border-blue-300 cursor-pointer hover:bg-blue-50 transition-colors">
                          <svg xmlns="http://www.w3.org/2000/svg" className="h-12 w-12 text-blue-400" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-8l-4-4m0 0L8 8m4-4v12" />
                          </svg>
                          <span className="mt-2 text-base font-medium">Upload ZIP Archive</span>
                          <span className="mt-1 text-xs text-gray-500">
                            Contains multiple X-ray images for batch analysis
                          </span>
                          <input
                            type="file"
                            accept=".zip"
                            onChange={handleBatchUpload}
                            className="hidden"
                          />
                        </label>
                      </div>

                      {batchFiles && (
                        <Card className="mb-6 border shadow-sm">
                          <CardHeader className="bg-gray-50 border-b">
                            <div className="flex justify-between items-center w-full">
                              <h3 className="text-md font-semibold">Batch Upload</h3>
                              <Chip size="sm" variant="flat">{batchFiles.length} file(s)</Chip>
                            </div>
                          </CardHeader>
                          <CardBody className="text-center py-6">
                            <svg xmlns="http://www.w3.org/2000/svg" className="h-12 w-12 mx-auto text-blue-500" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M5 8h14M5 8a2 2 0 110-4h14a2 2 0 110 4M5 8v10a2 2 0 002 2h10a2 2 0 002-2V8m-9 4h4" />
                            </svg>
                            <p className="mt-3 font-medium">Archive ready for processing</p>
                            <p className="text-sm text-gray-500 mt-1">
                              Archive name: {batchFiles[0]?.name}
                            </p>
                          </CardBody>
                          <CardFooter className="flex justify-center border-t bg-gray-50">
                            <Button
                              isDisabled={loading}
                              isLoading={loading}
                              onPress={predictBatchImages}
                              color="primary"
                              className="px-8 font-medium"
                              startContent={
                                <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5" viewBox="0 0 20 20" fill="currentColor">
                                  <path fillRule="evenodd" d="M7 2a1 1 0 00-.707 1.707L7 4.414v3.758a1 1 0 01-.293.707l-4 4C.817 14.769 2.156 18 4.828 18h10.343c2.673 0 4.012-3.231 2.122-5.121l-4-4A1 1 0 0113 8.172V4.414l.707-.707A1 1 0 0013 2H7zm2 6.172V4h2v4.172a3 3 0 00.879 2.12l1.168 1.168a4 4 0 00-2.929.923L9 13.5V12a1 1 0 00-1-1H5.415l-.879.879a1 1 0 001.415 1.415L7 12.245V13.5l-1.116-1.116a2.5 2.5 0 00-3.54 3.536l.708.708A4 4 0 006.586 18H4.828c-1.132 0-1.3-1.098-.708-1.708l4-4a3 3 0 00.879-2.12z" clipRule="evenodd" />
                                </svg>
                              }
                            >
                              Process Batch
                            </Button>
                          </CardFooter>
                        </Card>
                      )}

                      {batchResults.length > 0 && (
                        <Card className="mt-8 border shadow-md">
                          <CardHeader className="bg-gradient-to-r from-purple-800 to-purple-600 text-white">
                            <div className="flex justify-between items-center w-full">
                              <h3 className="text-xl font-bold">Batch Analysis Results</h3>
                              <Chip color="secondary" variant="dot">
                                {batchResults.length} images analyzed
                              </Chip>
                            </div>
                          </CardHeader>
                          <CardBody>
                            <div className="bg-purple-50 p-4 rounded-lg mb-6">
                              <div className="flex items-center">
                                <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5 text-purple-700 mr-2" viewBox="0 0 20 20" fill="currentColor">
                                  <path d="M2 11a1 1 0 011-1h2a1 1 0 011 1v5a1 1 0 01-1 1H3a1 1 0 01-1-1v-5zM8 7a1 1 0 011-1h2a1 1 0 011 1v9a1 1 0 01-1 1H9a1 1 0 01-1-1V7zM14 4a1 1 0 011-1h2a1 1 0 011 1v12a1 1 0 01-1 1h-2a1 1 0 01-1-1V4z" />
                                </svg>
                                <span className="font-medium text-purple-900">Batch Summary</span>
                              </div>
                              <div className="grid grid-cols-1 sm:grid-cols-3 gap-4 mt-3">
                                <div className="bg-white p-3 rounded-md shadow-sm">
                                  <div className="text-sm text-gray-500">Total Images</div>
                                  <div className="text-2xl font-bold text-purple-800">{batchResults.length}</div>
                                </div>
                                <div className="bg-white p-3 rounded-md shadow-sm">
                                  <div className="text-sm text-gray-500">Abnormal Findings</div>
                                  <div className="text-2xl font-bold text-red-600">
                                    {batchResults.filter(r =>
                                      r.Pathology !== 'No Finding' &&
                                      parseFloat(r.Probability.replace('%', '')) > 50
                                    ).length}
                                  </div>
                                </div>
                                <div className="bg-white p-3 rounded-md shadow-sm">
                                  <div className="text-sm text-gray-500">Normal</div>
                                  <div className="text-2xl font-bold text-green-600">
                                    {batchResults.filter(r =>
                                      r.Pathology === 'No Finding' ||
                                      parseFloat(r.Probability.replace('%', '')) <= 50
                                    ).length}
                                  </div>
                                </div>
                              </div>
                            </div>

                            <h4 className="text-md font-semibold mb-4">Detailed Batch Results</h4>
                            <div className="rounded-lg border overflow-hidden">
                              <Table
                                aria-label="Batch Prediction Results"
                                selectionMode="none"
                                isHeaderSticky
                                classNames={{
                                  base: "max-h-[480px]",
                                }}
                              >
                                <TableHeader>
                                  <TableColumn>File</TableColumn>
                                  <TableColumn>Primary Finding</TableColumn>
                                  <TableColumn>Confidence</TableColumn>
                                  <TableColumn>Severity</TableColumn>
                                  <TableColumn>Actions</TableColumn>
                                </TableHeader>
                                <TableBody emptyContent="No results available">
                                  {batchResults.map((result, idx) => {
                                    const probability = parseFloat(result.Probability.replace('%', ''));
                                    let severity = "Normal";
                                    let chipColor = "success";

                                    if (probability >= 80) {
                                      severity = "Critical";
                                      chipColor = "danger";
                                    } else if (probability >= 60) {
                                      severity = "Significant";
                                      chipColor = "warning";
                                    } else if (probability >= 40) {
                                      severity = "Moderate";
                                      chipColor = "primary";
                                    } else if (probability >= 20) {
                                      severity = "Mild";
                                      chipColor = "success";
                                    }

                                    return (
                                      <TableRow key={idx}>
                                        <TableCell>
                                          <div className="flex items-center">
                                            <svg xmlns="http://www.w3.org/2000/svg" className="h-4 w-4 text-gray-500 mr-2" viewBox="0 0 20 20" fill="currentColor">
                                              <path fillRule="evenodd" d="M4 4a2 2 0 012-2h4.586A2 2 0 0112 2.586L15.414 6A2 2 0 0116 7.414V16a2 2 0 01-2 2H6a2 2 0 01-2-2V4zm2 6a1 1 0 011-1h6a1 1 0 110 2H7a1 1 0 01-1-1zm1 3a1 1 0 100 2h6a1 1 0 100-2H7z" clipRule="evenodd" />
                                            </svg>
                                            <span className="text-sm">{result.FileName}</span>
                                          </div>
                                        </TableCell>
                                        <TableCell>{result.Pathology}</TableCell>
                                        <TableCell>
                                          <div className="flex items-center gap-2">
                                            <Progress
                                              aria-label="Loading..."
                                              size="sm"
                                              value={probability}
                                              color={chipColor}
                                              className="max-w-md"
                                            />
                                            <span className="text-sm font-medium">{result.Probability}</span>
                                          </div>
                                        </TableCell>
                                        <TableCell>
                                          <Chip size="sm" color={chipColor} variant="flat">{severity}</Chip>
                                        </TableCell>
                                        <TableCell>
                                          <div className="flex gap-2">
                                            <Button size="sm" variant="light" color="primary">
                                              View
                                            </Button>
                                            <Button size="sm" variant="light" color="success">
                                              Report
                                            </Button>
                                          </div>
                                        </TableCell>
                                      </TableRow>
                                    );
                                  })}
                                </TableBody>
                              </Table>
                            </div>
                          </CardBody>
                          <CardFooter className="flex justify-center border-t pt-4 gap-2">
                            <Button
                              color="default"
                              variant="flat"
                              startContent={
                                <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5" viewBox="0 0 20 20" fill="currentColor">
                                  <path fillRule="evenodd" d="M3 17a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1zm3.293-7.707a1 1 0 011.414 0L9 10.586V3a1 1 0 112 0v7.586l1.293-1.293a1 1 0 111.414 1.414l-3 3a1 1 0 01-1.414 0l-3-3a1 1 0 010-1.414z" clipRule="evenodd" />
                                </svg>
                              }
                            >
                              Export All Results
                            </Button>
                            <Button
                              color="primary"
                              variant="flat"
                              startContent={
                                <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5" viewBox="0 0 20 20" fill="currentColor">
                                  <path d="M4 3a2 2 0 100 4h12a2 2 0 100-4H4z" />
                                  <path fillRule="evenodd" d="M3 8h14v7a2 2 0 01-2 2H5a2 2 0 01-2-2V8zm5 3a1 1 0 011-1h2a1 1 0 110 2H9a1 1 0 01-1-1z" clipRule="evenodd" />
                                </svg>
                              }
                            >
                              Print Summary
                            </Button>
                          </CardFooter>
                        </Card>
                      )}
                    </div>
                  </div>
                </Tab>
              </Tabs>
            </div>
          </CardBody>
        </Card>

        {error && (
          <div className="bg-red-50 border border-red-200 text-red-700 px-4 py-3 rounded-lg flex items-center mb-6">
            <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5 mr-2 text-red-500" viewBox="0 0 20 20" fill="currentColor">
              <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7 4a1 1 0 11-2 0 1 1 0 012 0zm-1-9a1 1 0 00-1 1v4a1 1 0 102 0V6a1 1 0 00-1-1z" clipRule="evenodd" />
            </svg>
            {error}
          </div>
        )}

        <div className="grid grid-cols-1 md:grid-cols-3 gap-6 mb-8">
          <Card className="border shadow-sm">
            <CardBody className="flex items-start">
              <div className="bg-blue-100 p-3 rounded-full mr-4">
                <svg xmlns="http://www.w3.org/2000/svg" className="h-6 w-6 text-blue-700" viewBox="0 0 20 20" fill="currentColor">
                  <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2h-1V9a1 1 0 00-1-1z" clipRule="evenodd" />
                </svg>
              </div>
              <div>
                <h3 className="font-semibold text-lg mb-1">Accuracy Information</h3>
                <p className="text-sm text-gray-600">
                  The AI model has a 92.8% accuracy rate on validated datasets. Results should be confirmed by a medical professional.
                </p>
              </div>
            </CardBody>
          </Card>

          <Card className="border shadow-sm">
            <CardBody className="flex items-start">
              <div className="bg-purple-100 p-3 rounded-full mr-4">
                <svg xmlns="http://www.w3.org/2000/svg" className="h-6 w-6 text-purple-700" viewBox="0 0 20 20" fill="currentColor">
                  <path d="M9 4.804A7.968 7.968 0 005.5 4c-1.255 0-2.443.29-3.5.804v10A7.969 7.969 0 015.5 14c1.669 0 3.218.51 4.5 1.385A7.962 7.962 0 0114.5 14c1.255 0 2.443.29 3.5.804v-10A7.968 7.968 0 0014.5 4c-1.255 0-2.443.29-3.5.804V12a1 1 0 11-2 0V4.804z" />
                </svg>
              </div>
              <div>
                <h3 className="font-semibold text-lg mb-1">Research Reference</h3>
                <p className="text-sm text-gray-600">
                  Based on deep learning models trained on over 200,000 annotated radiographic images from clinical studies.
                </p>
              </div>
            </CardBody>
          </Card>

          <Card className="border shadow-sm">
            <CardBody className="flex items-start">
              <div className="bg-green-100 p-3 rounded-full mr-4">
                <svg xmlns="http://www.w3.org/2000/svg" className="h-6 w-6 text-green-700" viewBox="0 0 20 20" fill="currentColor">
                  <path fillRule="evenodd" d="M2.166 4.999A11.954 11.954 0 0010 1.944 11.954 11.954 0 0017.834 5c.11.65.166 1.32.166 2.001 0 5.225-3.34 9.67-8 11.317C5.34 16.67 2 12.225 2 7c0-.682.057-1.35.166-2.001zm11.541 3.708a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z" clipRule="evenodd" />
                </svg>
              </div>
              <div>
                <h3 className="font-semibold text-lg mb-1">Data Privacy</h3>
                <p className="text-sm text-gray-600">
                  All uploaded images are processed securely and confidentially in compliance with HIPAA regulations.
                </p>
              </div>
            </CardBody>
          </Card>
        </div>
      </main>

      <footer className="bg-blue-900 text-white py-6 px-8 mt-auto">
        <div className="max-w-7xl mx-auto">
          <div className="flex flex-col md:flex-row justify-between items-center">
            <div className="mb-4 md:mb-0">
              <div className="flex items-center">
                <svg xmlns="http://www.w3.org/2000/svg" className="h-6 w-6 mr-2" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 3v2m6-2v2M9 19v2m6-2v2M5 9H3m2 6H3m18-6h-2m2 6h-2M7 19h10a2 2 0 002-2V7a2 2 0 00-2-2H7a2 2 0 00-2 2v10a2 2 0 002 2zM9 9h6v6H9V9z" />
                </svg>
                <h2 className="text-lg font-bold">Medical Imaging AI Platform</h2>
              </div>
              <p className="text-blue-200 text-sm mt-1">Advanced diagnostics powered by AI</p>
            </div>
            <div className="flex flex-col md:flex-row md:items-center text-center md:text-left">
              <div className="text-sm text-blue-200 mr-6 mb-2 md:mb-0">
                <span>Model: ChestX-Net v3.2</span>
              </div>
              <div className="text-sm text-blue-200 mr-6 mb-2 md:mb-0">
                <span>Last updated: April 2025</span>
              </div>
              <Button size="sm" color="primary" variant="solid" className="bg-blue-700">
                Documentation
              </Button>
            </div>
          </div>
          <div className="mt-6 pt-6 border-t border-blue-800 text-center md:text-left text-xs text-blue-300">
            <p>© 2025 Medical Imaging AI. For clinical decision support only. Not intended to replace professional medical diagnosis.</p>
          </div>
        </div>
      </footer>
    </div>
  );
}