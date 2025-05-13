'use client';

import React, { useState } from 'react';
import { Button } from '@nextui-org/button';
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer
} from 'recharts';
import { AlertCircle, UploadCloud, FileImage, Activity } from 'lucide-react';

export default function SkinCancerAnalysis() {
  const [uploadedFile, setUploadedFile] = useState(null);
  const [selectedFile, setSelectedFile] = useState(null);
  const [prediction, setPrediction] = useState(null);
  const [probabilities, setProbabilities] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [showDemoData, setShowDemoData] = useState(false);

  const handleFileUpload = (event) => {
    const file = event.target.files?.[0];
    if (!file) return;

    setError(null);
    setUploadedFile(file);
    setSelectedFile(file.name);
    setPrediction(null);
    setProbabilities(null);
    setShowDemoData(false);
  };

  const handlePrediction = async () => {
    if (!uploadedFile) {
      setError('Please upload a dermatoscopic image.');
      return;
    }

    try {
      setLoading(true);

      // In a real application, this would send data to the server
      // For demonstration purposes only, we'll simulate a response
      setTimeout(() => {
        setLoading(false);
        setShowDemoData(true);
      }, 2000);

    } catch (err) {
      setError('Error occurred during image analysis.');
      setLoading(false);
    }
  };

  // Demo data - would be replaced by actual API response
  const demoData = {
    prediction: "Benign nevus",
    probabilities: {
      "Benign nevus": 0.82,
      "Melanoma": 0.09,
      "Seborrheic keratosis": 0.05,
      "Basal cell carcinoma": 0.03,
      "Actinic keratosis": 0.01
    }
  };

  const getChartData = () => {
    if (showDemoData) {
      return Object.entries(demoData.probabilities)
        .map(([key, value]) => ({
          name: key,
          probability: (value * 100).toFixed(1)
        }))
        .sort((a, b) => b.probability - a.probability);
    }
    return [];
  };

  // Color mapping for different diagnosis types
  const getBarColor = (entry) => {
    const name = entry.name;
    if (name.includes("Melanoma")) return "#e74c3c";
    if (name.includes("Basal cell") || name.includes("Squamous cell")) return "#e67e22";
    if (name.includes("Actinic") || name.includes("Seborrheic")) return "#f39c12";
    return "#2ecc71"; // Benign by default
  };

  return (
    <div className="bg-gray-50 min-h-screen">
      <header className="bg-white shadow-md">
        <div className="max-w-7xl mx-auto py-4 px-6">
          <div className="flex items-center justify-between">
            <div className="flex items-center space-x-2">
              <Activity size={28} className="text-blue-600" />
              <h1 className="text-2xl font-bold text-gray-800">Dermatology  Analysis System</h1>
            </div>
            <div className="text-sm text-gray-500">v1.0.0</div>
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto py-8 px-6">
        <div className="bg-white shadow-lg rounded-lg overflow-hidden mb-8">
          <div className="bg-blue-600 py-4 px-6">
            <h2 className="text-white text-lg font-semibold flex items-center">
              <FileImage className="mr-2" size={20} />
              Image Analysis Module
            </h2>
          </div>

          <div className="p-6">
            <div className="mb-8">
              <h3 className="text-gray-700 font-medium mb-4">Upload Dermatoscopic Image</h3>
              <div className="border-2 border-dashed border-gray-300 rounded-lg p-8 text-center hover:border-blue-500 transition-colors">
                <UploadCloud className="mx-auto h-12 w-12 text-gray-400" />
                <p className="mt-2 text-sm text-gray-500">Supported formats: JPG, PNG (up to 10MB)</p>
                <input
                  type="file"
                  accept="image/*"
                  onChange={handleFileUpload}
                  className="hidden"
                  id="file-upload"
                />
                <label
                  htmlFor="file-upload"
                  className="mt-4 inline-flex items-center px-4 py-2 border border-transparent text-sm font-medium rounded-md shadow-sm text-white bg-blue-600 hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-blue-500 cursor-pointer"
                >
                  Select Image
                </label>
              </div>

              {selectedFile && (
                <div className="mt-4 bg-blue-50 p-4 rounded-md flex items-center">
                  <FileImage size={20} className="text-blue-600 mr-2" />
                  <span className="text-sm">{selectedFile}</span>
                </div>
              )}

              <div className="mt-6">
                <Button
                  isDisabled={loading || !uploadedFile}
                  isLoading={loading}
                  onPress={handlePrediction}
                  className="bg-blue-600 hover:bg-blue-700 text-white font-medium py-2 px-6 rounded-md"
                >
                  {loading ? "Analyzing..." : "Analyze Image"}
                </Button>
              </div>

              {error && (
                <div className="mt-4 bg-red-50 p-4 rounded-md flex items-center">
                  <AlertCircle size={20} className="text-red-600 mr-2" />
                  <span className="text-sm text-red-600">{error}</span>
                </div>
              )}
            </div>

            {showDemoData && (
              <div className="mt-8 border-t pt-8">
                <div className="grid md:grid-cols-2 gap-8">
                  <div>
                    <h3 className="text-lg font-semibold mb-4 text-gray-800">Analysis Results</h3>
                    <div className="bg-gray-50 p-6 rounded-lg">
                      <div className="mb-6">
                        <h4 className="text-sm font-medium text-gray-500 mb-1">Primary Diagnosis</h4>
                        <p className="text-2xl font-bold text-green-600">{demoData.prediction}</p>
                      </div>

                      <h4 className="text-sm font-medium text-gray-500 mb-3">Differential Diagnosis</h4>
                      <div className="space-y-3">
                        {Object.entries(demoData.probabilities)
                          .sort(([, a], [, b]) => b - a)
                          .map(([diagnosis, probability]) => (
                            <div key={diagnosis} className="flex items-center justify-between">
                              <span className="text-sm font-medium">{diagnosis}</span>
                              <div className="flex items-center">
                                <div className="w-36 bg-gray-200 rounded-full h-2 mr-2">
                                  <div
                                    className={`h-2 rounded-full ${
                                      diagnosis === "Melanoma" ? "bg-red-500" : 
                                      diagnosis === "Basal cell carcinoma" ? "bg-orange-500" :
                                      diagnosis === "Actinic keratosis" || diagnosis === "Seborrheic keratosis" ? "bg-yellow-500" :
                                      "bg-green-500"
                                    }`}
                                    style={{ width: `${(probability * 100)}%` }}
                                  ></div>
                                </div>
                                <span className="text-sm font-medium">{(probability * 100).toFixed(1)}%</span>
                              </div>
                            </div>
                          ))
                        }
                      </div>

                      <div className="mt-6 text-xs text-gray-500">
                        <p>Note: This is a clinical decision support tool. Results should be interpreted by a qualified healthcare professional.</p>
                      </div>
                    </div>
                  </div>

                  <div>
                    <h3 className="text-lg font-semibold mb-4 text-gray-800">Probability Distribution</h3>
                    <div className="bg-gray-50 p-4 rounded-lg" style={{ height: '350px' }}>
                      <ResponsiveContainer width="100%" height="100%">
                        <BarChart
                          data={getChartData()}
                          layout="vertical"
                          margin={{ top: 5, right: 30, left: 20, bottom: 5 }}
                        >
                          <CartesianGrid strokeDasharray="3 3" horizontal={false} />
                          <XAxis type="number" domain={[0, 100]} unit="%" />
                          <YAxis dataKey="name" type="category" width={150} />
                          <Tooltip
                            formatter={(value) => [`${value}%`, 'Probability']}
                            contentStyle={{ borderRadius: '4px' }}
                          />
                          <Bar
                            dataKey="probability"
                            name="Probability"
                            barSize={24}
                            fill="#4f46e5"
                            radius={[0, 4, 4, 0]}
                          />
                        </BarChart>
                      </ResponsiveContainer>
                    </div>
                  </div>
                </div>
              </div>
            )}
          </div>
        </div>

        <div className="bg-white shadow-lg rounded-lg overflow-hidden mb-8">
          <div className="bg-blue-600 py-4 px-6">
            <h2 className="text-white text-lg font-semibold">About This Tool</h2>
          </div>
          <div className="p-6">
            <p className="text-gray-600 text-sm mb-4">
              This clinical decision support tool uses deep learning algorithms to analyze dermatoscopic images and provide probabilistic assessments of skin lesions. The system is trained on a diverse dataset of clinically validated images.
            </p>
            <p className="text-gray-600 text-sm">
              <strong>Important:</strong> This tool is designed to assist healthcare professionals and should not replace clinical judgment. All results should be interpreted within the appropriate clinical context.
            </p>
          </div>
        </div>
      </main>

      <footer className="bg-gray-800 text-gray-300 py-6">
        <div className="max-w-7xl mx-auto px-6">
          <div className="flex flex-col md:flex-row justify-between items-center">
            <div className="mb-4 md:mb-0">
              <div className="flex items-center space-x-2">
                <Activity size={20} />
                <span className="font-semibold">Dermatology  Analysis System</span>
              </div>
              <p className="text-xs mt-1">For research and educational purposes only</p>
            </div>
            <div className="text-xs">
              © {new Date().getFullYear()} Medical Imaging Research Group
            </div>
          </div>
        </div>
      </footer>
    </div>
  );
}