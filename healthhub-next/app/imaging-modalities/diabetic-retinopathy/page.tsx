"use client";

import React, { useState } from "react";
import axios from "axios";
import Image from "next/image";
import { Card, CardBody, CardFooter } from "@nextui-org/card";
import { Button } from "@nextui-org/button";
import { Progress } from "@nextui-org/progress";
import { Tabs, Tab } from "@nextui-org/tabs";
import { Tooltip } from "@nextui-org/tooltip";

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
  const [activeTab, setActiveTab] = useState<string>("single");
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
    formData.append("file", file);

    try {
      setLoading(true);
      const response = await axios.post<PredictionResult>(
        "http://localhost:8000/imaging/diabetic-retinopathy/predict",
        formData,
        {
          headers: { "Content-Type": "multipart/form-data" },
        }
      );
      setPredictionResult(response.data);
      console.log(response.data);
    } catch (err) {
      setError("An error occurred while processing the image. Please try again.");
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
    files.forEach((file) => formData.append("files", file));

    try {
      setBatchLoading(true);
      const response = await axios.post<{ predictions: BatchPrediction[] }>(
        "http://localhost:8000/imaging/diabetic-retinopathy/batch-predict",
        formData,
        {
          headers: { "Content-Type": "multipart/form-data" },
        }
      );
      setBatchResults(response.data.predictions);
    } catch (err) {
      setError("An error occurred while processing the batch upload. Please try again.");
    } finally {
      setBatchLoading(false);
    }
  };

  return (
    <div className="max-w-7xl mx-auto p-8">
      <h1 className="text-3xl font-bold text-center mb-4">Diabetic Retinopathy Classification</h1>
      <p className="text-lg text-gray-600 text-center mb-8">
        Upload your retinal images to classify them for diabetic retinopathy using advanced deep learning models.
      </p>

      <Tabs activeTab={activeTab} onChange={setActiveTab} className="mb-8">
        {/* Single Image Tab */}
        <Tab title="Single Image" key="single">
          <div>
            <Card className="mb-6">
              <CardBody>
                <h2 className="text-xl font-bold mb-4">Upload a Single Image</h2>
                <input
                  type="file"
                  accept="image/*"
                  onChange={handleImageUpload}
                  className="block w-full max-w-sm text-sm text-gray-500 file:mr-4 file:py-2 file:px-4 file:rounded-full file:border-0 file:bg-indigo-50 file:text-indigo-700 hover:file:bg-indigo-100"
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
                    className="object-cover rounded-lg"
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
                    {loading ? "Processing..." : "Predict"}
                  </Button>
                </CardFooter>
              </Card>
            )}

            {predictionResult && (
              <Card className="mb-6">
                <CardBody>
                  <h2 className="text-xl font-bold text-center mb-2">Prediction Results</h2>
                  <p
                    className={`text-lg font-semibold ${
                      predictionResult.result === "Severe" ? "text-red-600" : "text-green-600"
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
        </Tab>

        {/* Batch Upload Tab */}
        <Tab title="Batch Upload" key="batch">
          <div>
            <Card className="mb-6">
              <CardBody>
                <h2 className="text-xl font-bold mb-4">Batch Upload (ZIP File)</h2>
                <input
                  type="file"
                  accept=".zip"
                  onChange={handleBatchUpload}
                  className="block w-full max-w-sm text-sm text-gray-500 file:mr-4 file:py-2 file:px-4 file:rounded-full file:border-0 file:bg-indigo-50 file:text-indigo-700 hover:file:bg-indigo-100"
                />
              </CardBody>
            </Card>

            {batchLoading && <Progress value={100} label="Processing Batch..." />}

            {batchResults.length > 0 && (
              <Card className="mb-6">
                <CardBody>
                  <h2 className="text-xl font-bold text-center mb-4">Batch Prediction Results</h2>
                  <ul>
                    {batchResults.map((result, index) => (
                      <li key={index} className="text-lg mb-2">
                        <b>Image {index + 1}:</b> {result.result} ({(result.confidence * 100).toFixed(2)}%)
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
        </Tab>
      </Tabs>

      <Tooltip content="Sample Images Coming Soon" placement="bottom">
        <Button isDisabled className="mt-6 mx-auto block">
          View Sample Images
        </Button>
      </Tooltip>
    </div>
  );
}
