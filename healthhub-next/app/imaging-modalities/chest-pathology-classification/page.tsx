'use client';

import React, { useState } from 'react';
import axios from 'axios';
import { Card, CardBody, CardFooter } from '@nextui-org/card';
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
} from 'recharts';
import Image from 'next/image';
import { Button } from '@nextui-org/button';
import { Switch } from '@nextui-org/switch';
import {cn} from "@nextui-org/react"
// import cn from 'classnames';

export default function ChestXRayClassification() {
  const [uploadedImage, setUploadedImage] = useState<string | null>(null);
  const [file, setFile] = useState<File | null>(null);
  const [isDicomFile, setIsDicomFile] = useState<boolean>(false);
  const [predictionResults, setPredictionResults] = useState([]);
  const [batchFiles, setBatchFiles] = useState<FileList | null>(null);
  const [batchResults, setBatchResults] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [activeTab, setActiveTab] = useState<string>('single');

  const handleSingleImageUpload = (
    event: React.ChangeEvent<HTMLInputElement>
  ) => {
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
    return predictionResults.map((result: any) => ({
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

  return (
    <div className="max-w-7xl mx-auto p-8">
      <h1 className="text-4xl font-extrabold text-center mb-4">
        🩺 Chest X-Ray Classification
      </h1>
      <p className="text-lg text-gray-600 text-center mb-6">
        Upload chest X-ray images (JPEG, PNG, or DICOM) to identify potential
        pathologies using our AI-powered model.
      </p>

      {/* Custom Switch to toggle between Single Image and Batch Prediction */}
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
              {activeTab === 'batch' ? 'Batch Prediction' : 'Single Image'}
            </p>
            <p className="text-tiny text-default-400">
              {activeTab === 'batch'
                ? 'Upload a ZIP file containing multiple chest X-ray images.'
                : 'Upload a single chest X-ray image (JPEG, PNG, or DICOM).'}
            </p>
          </div>
        </Switch>
      </div>

      {/* Section Heading */}
      <h2 className="text-2xl font-bold mb-6 text-center">
        {activeTab === 'single' ? 'Upload a Chest X-Ray' : 'Batch Prediction'}
      </h2>

      {activeTab === 'single' && (
        <section>
          <p className="text-center mb-4 text-gray-500">
            Supported formats: JPEG, PNG, or DICOM (.dcm)
          </p>
          <div className="flex justify-center">
            <input
              type="file"
              accept="image/*, .dcm"
              onChange={handleSingleImageUpload}
              className="block w-full max-w-md text-sm text-gray-500
                  file:mr-4 file:py-2 file:px-4 file:rounded-full
                  file:border-0 file:bg-indigo-600 file:text-white
                  hover:file:bg-indigo-700 mb-4 cursor-pointer"
            />
          </div>
          {(uploadedImage || file) && (
            <Card className="mb-6">
              <CardBody>
                {isDicomFile ? (
                  <div className="text-center py-20">
                    <p className="text-lg font-semibold">
                      DICOM file "{file?.name}" has been uploaded successfully.
                    </p>
                    <p className="text-gray-500 mt-2">
                      DICOM images cannot be previewed, but you can proceed to
                      prediction.
                    </p>
                  </div>
                ) : (
                  <div className="flex justify-center">
                    <Image
                      src={uploadedImage!}
                      alt="Uploaded X-Ray"
                      width={600}
                      height={400}
                      className="rounded-lg"
                    />
                  </div>
                )}
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
              <h3 className="text-2xl font-bold mb-4 text-center">
                Prediction Results
              </h3>
              <div className="grid grid-cols-1 lg:grid-cols-2 gap-8 mb-6">
                {/* Bar Chart */}
                <div className="h-96">
                  <ResponsiveContainer>
                    <BarChart data={transformPredictionDataForChart()}>
                      <CartesianGrid strokeDasharray="3 3" />
                      <XAxis dataKey="name" />
                      <YAxis />
                      <Tooltip />
                      <Bar
                        dataKey="probability"
                        fill="#8884d8"
                        label={{ position: 'top' }}
                      >
                        {transformPredictionDataForChart().map(
                          (entry, index) => (
                            <Cell
                              key={`cell-${index}`}
                              fill={
                                topPrediction &&
                                entry.name === topPrediction.Pathology
                                  ? '#FF0000'
                                  : COLORS[index % COLORS.length]
                              }
                            />
                          )
                        )}
                      </Bar>
                    </BarChart>
                  </ResponsiveContainer>
                </div>

                {/* Donut Chart */}
                <div className="h-96">
                  <ResponsiveContainer>
                    <PieChart>
                      <Pie
                        data={transformPredictionDataForChart()}
                        dataKey="probability"
                        nameKey="name"
                        cx="50%"
                        cy="50%"
                        outerRadius={120}
                        innerRadius={60}
                        fill="#8884d8"
                        label={({ name, probability }) =>
                          `${name}: ${probability.toFixed(1)}%`
                        }
                        labelLine={false}
                      >
                        {transformPredictionDataForChart().map(
                          (entry, index) => (
                            <Cell
                              key={`cell-${index}`}
                              fill={
                                topPrediction &&
                                entry.name === topPrediction.Pathology
                                  ? '#FF0000'
                                  : COLORS[index % COLORS.length]
                              }
                            />
                          )
                        )}
                      </Pie>
                      <Tooltip />
                      <Legend verticalAlign="bottom" height={36} />
                    </PieChart>
                  </ResponsiveContainer>
                </div>
              </div>

              <Table
                aria-label="Prediction Results"
                className="mt-6"
                striped
                sticked
              >
                <TableHeader>
                  <TableColumn>Pathology</TableColumn>
                  <TableColumn>Probability</TableColumn>
                </TableHeader>
                <TableBody>
                  {predictionResults.map((result: any, idx: number) => (
                    <TableRow
                      key={idx}
                      className={
                        topPrediction &&
                        result.Pathology === topPrediction.Pathology
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

          {error && (
            <p className="text-red-600 mt-4 text-center font-semibold">
              {error}
            </p>
          )}
        </section>
      )}

      {activeTab === 'batch' && (
        <section>
          <p className="text-center mb-4 text-gray-500">
            Upload a ZIP file containing multiple chest X-ray images.
          </p>
          <div className="flex justify-center">
            <input
              type="file"
              accept=".zip"
              onChange={handleBatchUpload}
              className="block w-full max-w-md text-sm text-gray-500
                  file:mr-4 file:py-2 file:px-4 file:rounded-full
                  file:border-0 file:bg-green-600 file:text-white
                  hover:file:bg-green-700 mb-4 cursor-pointer"
            />
          </div>
          <div className="flex justify-center">
            <Button
              isDisabled={loading || !batchFiles}
              isLoading={loading}
              onPress={predictBatchImages}
              color="primary"
              className="mt-4"
            >
              Predict Batch
            </Button>
          </div>

          {batchResults.length > 0 && (
            <div className="mt-8">
              <h3 className="text-2xl font-bold mb-4 text-center">
                Batch Prediction Results
              </h3>
              <Table
                aria-label="Batch Prediction Results"
                striped
                sticked
                className="mt-6"
              >
                <TableHeader>
                  <TableColumn>File Name</TableColumn>
                  <TableColumn>Pathology</TableColumn>
                  <TableColumn>Probability</TableColumn>
                </TableHeader>
                <TableBody>
                  {batchResults.map((result: any, idx: number) => (
                    <TableRow key={idx}>
                      <TableCell>{result.FileName}</TableCell>
                      <TableCell>{result.Pathology}</TableCell>
                      <TableCell>{result.Probability}</TableCell>
                    </TableRow>
                  ))}
                </TableBody>
              </Table>
            </div>
          )}

          {error && (
            <p className="text-red-600 mt-4 text-center font-semibold">
              {error}
            </p>
          )}
        </section>
      )}
    </div>
  );
}
