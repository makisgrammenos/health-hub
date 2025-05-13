'use client'
import React from "react";
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from "recharts";
import { Brain, Heart, Activity, FileBarChart, PieChart, Zap, Database, GitMerge } from "lucide-react";

// Mock data for charts
const analyticsData = [
  { name: "Jan", cases: 400, predictions: 240 },
  { name: "Feb", cases: 300, predictions: 250 },
  { name: "Mar", cases: 520, predictions: 490 },
  { name: "Apr", cases: 480, predictions: 400 },
  { name: "May", cases: 390, predictions: 350 },
  { name: "Jun", cases: 600, predictions: 570 }
];

const ServiceCard = ({ service }) => {
  return (
    <a
      href={service.path}
      className="bg-gray-50 border border-gray-200 rounded-lg p-4 hover:bg-blue-50 transition-colors cursor-pointer flex flex-col h-full"
    >
      <div className="h-40 bg-gray-100 rounded-md flex items-center justify-center mb-4 overflow-hidden">
        {getServiceIcon(service.title)}
      </div>
      <h3 className="text-lg font-semibold text-blue-800 mb-2">{service.title}</h3>
      <p className="text-gray-600 text-sm flex-grow">{service.description}</p>
      <div className="mt-4 flex justify-between items-center">
        <span className="text-xs text-gray-500">ML-powered analysis</span>
        <span className="bg-blue-100 text-blue-800 text-xs px-2 py-1 rounded-full">Clinical</span>
      </div>
    </a>
  );
};

const getServiceIcon = (title) => {
  if (title.includes("Heart")) return <Heart size={64} className="text-red-500" />;
  if (title.includes("Plot")) return <FileBarChart size={64} className="text-blue-500" />;
  if (title.includes("Dimensionality")) return <GitMerge size={64} className="text-purple-500" />;
  if (title.includes("Classification")) return <Brain size={64} className="text-green-500" />;
  if (title.includes("Clustering")) return <Database size={64} className="text-yellow-500" />;
  if (title.includes("Feature")) return <PieChart size={64} className="text-indigo-500" />;
  return <Activity size={64} className="text-gray-500" />;
};

export default function ScientificInterface() {
  const services = [
    {
      title: "Heart Data preprocessing",
      description: "Leverage AI to detect diabetic retinopathy with precision.",
      path: "/clinical-tools/basic-processing",
    },
    {
      title: "Plot heart disease data",
      description: "Create various plots to visualize heart disease data.",
      path: "/clinical-tools/plots",
    },
    {
      title: "Dimensionality Reduction",
      description: "Reduce the number of features in a dataset while preserving its structure using various techniques.",
      path: "/clinical-tools/dimensionality-reduction",
    },
    {
      title: "Classification using Machine Learning models",
      description: "Train and evaluate machine learning models to classify heart disease data.",
      path: "/clinical-tools/classification",
    },
    {
      title: "Clustering",
      description: "Group similar data points together using clustering algorithms.",
      path: "/clinical-tools/clustering",
    },
    {
      title: "Feature Selection and Importance",
      description: "Identify the most important features in a dataset using feature selection techniques.",
      path: "/clinical-tools/feature-importance",
    },
  ];

  return (
    <div className="min-h-screen bg-gray-50">
      {/* Header Bar */}
      <header className="bg-blue-900 text-white py-3 px-6 flex items-center justify-between shadow-md">
        <div className="flex items-center">
          <Activity className="mr-2" size={24} />
          <h1 className="text-xl font-mono font-bold">CardioML Research Suite</h1>
        </div>
        <div className="flex items-center space-x-2 text-sm">
          <span className="px-3 py-1 bg-blue-800 rounded-md">Data Import</span>
          <span className="px-3 py-1 bg-blue-800 rounded-md">Analysis</span>
          <span className="px-3 py-1 bg-blue-800 rounded-md">Publication</span>
          <span className="px-3 py-1 bg-green-600 rounded-md flex items-center">
            <Zap size={16} className="mr-1" />
            Online
          </span>
        </div>
      </header>

      <div className="container mx-auto px-4 py-6">
        {/* Dashboard Header */}
        <div className="mb-8">
          <h1 className="text-3xl font-bold text-blue-900 mb-2">Clinical Data Analysis Platform</h1>
          <p className="text-gray-600">Advanced machine learning tools for cardiac research and diagnostics</p>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mt-6">
            <div className="bg-gradient-to-r from-blue-600 to-blue-700 rounded-lg p-4 text-white flex items-center shadow-md">
              <div className="bg-white/20 p-3 rounded-lg mr-4">
                <Database size={24} />
              </div>
              <div>
                <div className="text-sm opacity-80">Imported Datasets</div>
                <div className="text-2xl font-semibold">10</div>
              </div>
            </div>

            <div className="bg-gradient-to-r from-green-600 to-green-700 rounded-lg p-4 text-white flex items-center shadow-md">
              <div className="bg-white/20 p-3 rounded-lg mr-4">
                <Heart size={24} />
              </div>
              <div>
                <div className="text-sm opacity-80">Cases Processed</div>
                <div className="text-2xl font-semibold">127</div>
              </div>
            </div>

            <div className="bg-gradient-to-r from-purple-600 to-purple-700 rounded-lg p-4 text-white flex items-center shadow-md">
              <div className="bg-white/20 p-3 rounded-lg mr-4">
                <Brain size={24} />
              </div>
              <div>
                <div className="text-sm opacity-80">Average AUC</div>
                <div className="text-2xl font-semibold">94.4%</div>
              </div>
            </div>
          </div>
        </div>

        {/* Tool Selection */}
        <div className="bg-white rounded-lg border border-gray-200 p-6 shadow-sm mb-8">
          <div className="flex justify-between items-center mb-6">
            <h2 className="text-xl font-bold text-gray-800">Cardiac ML Tools</h2>
            <div className="text-sm text-blue-600 font-medium">System v2.4.1</div>
          </div>

          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
            {services.map((service, index) => (
              <ServiceCard
                key={index}
                service={service}
              />
            ))}
          </div>
        </div>

        {/* Recent Activity */}
        {/*<div className="bg-white rounded-lg border border-gray-200 p-6 shadow-sm">*/}
        {/*  <h2 className="text-xl font-bold text-gray-800 mb-4">Research Metrics</h2>*/}
        {/*  <ResponsiveContainer width="100%" height={300}>*/}
        {/*    <BarChart data={analyticsData}>*/}
        {/*      <CartesianGrid strokeDasharray="3 3" />*/}
        {/*      <XAxis dataKey="name" />*/}
        {/*      <YAxis />*/}
        {/*      <Tooltip />*/}
        {/*      <Legend />*/}
        {/*      <Bar dataKey="cases" fill="#3182CE" name="Clinical Cases" />*/}
        {/*      <Bar dataKey="predictions" fill="#48BB78" name="ML Predictions" />*/}
        {/*    </BarChart>*/}
        {/*  </ResponsiveContainer>*/}
        {/*</div>*/}
      </div>

      {/* Footer */}
      <footer className="bg-gray-100 border-t border-gray-200 py-4 px-6 mt-8">
        <div className="container mx-auto flex flex-col md:flex-row justify-between items-center text-sm text-gray-600">
          <div>CardioML Research Suite © 2025 | Data Analytics Platform v2.4.1</div>
          <div className="flex space-x-4 mt-2 md:mt-0">
            <span>Documentation</span>
            <span>Changelog</span>
            <span>Research Papers</span>
          </div>
        </div>
      </footer>
    </div>
  );
}