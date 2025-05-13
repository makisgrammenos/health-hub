'use client';

import { useState, useEffect } from "react";
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';
import { ArrowUpRight, AlertTriangle, CheckCircle, Database, BarChart2, FileText } from "lucide-react";

const ClassificationDemo = () => {
  const [metrics, setMetrics] = useState([]);
  const [predictions, setPredictions] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [activeTab, setActiveTab] = useState('metrics');
  const [visibleColumns, setVisibleColumns] = useState({});

  useEffect(() => {
    const fetchData = async () => {
      try {
        const response = await fetch("http://localhost:8000/clinical/classification");
        if (!response.ok) throw new Error("Failed to fetch classification data");
        const data = await response.json();
        setMetrics(data.metrics);
        setPredictions(data.predictions);

        // Set initial visible columns
        if (data.metrics.length > 0) {
          const metricColumns = Object.keys(data.metrics[0]);
          const initialVisibleColumns = {};
          metricColumns.forEach(col => {
            initialVisibleColumns[col] = true;
          });
          setVisibleColumns(initialVisibleColumns);
        }
      } catch (error) {
        console.error("Error fetching data:", error);
        setError(error.message);
      } finally {
        setLoading(false);
      }
    };

    fetchData();
  }, []);

  const toggleColumnVisibility = (column) => {
    setVisibleColumns(prev => ({
      ...prev,
      [column]: !prev[column]
    }));
  };

  // Prepare chart data - find appropriate keys dynamically
  const prepareChartData = () => {
    if (!metrics.length) return [];

    // Find likely column names
    const keys = Object.keys(metrics[0]);
    const modelNameKey = keys.find(key =>
      key.toLowerCase().includes('model') || key.toLowerCase().includes('name')
    ) || keys[0];

    const accuracyKey = keys.find(key => key.toLowerCase().includes('accuracy'));
    const precisionKey = keys.find(key => key.toLowerCase().includes('precision'));
    const recallKey = keys.find(key => key.toLowerCase().includes('recall'));
    const f1Key = keys.find(key =>
      key.toLowerCase().includes('f1') || key.toLowerCase().includes('f_1') || key.toLowerCase().includes('f-score')
    );

    return metrics.map(metric => {
      const chartPoint = {
        name: String(metric[modelNameKey] || "Model"),
      };

      if (accuracyKey) chartPoint.accuracy = Number(metric[accuracyKey] || 0);
      if (precisionKey) chartPoint.precision = Number(metric[precisionKey] || 0);
      if (recallKey) chartPoint.recall = Number(metric[recallKey] || 0);
      if (f1Key) chartPoint.f1 = Number(metric[f1Key] || 0);

      return chartPoint;
    });
  };

  const chartData = prepareChartData();

  // Calculate summary metrics
  const calculateSummary = () => {
    if (!metrics.length) return {
      bestModel: "No data",
      bestAccuracy: 0,
      avgAccuracy: 0,
      totalModels: 0
    };

    // Define expected metric column names and fallbacks
    const modelNameKey = Object.keys(metrics[0]).find(key =>
      key.toLowerCase().includes('model') || key.toLowerCase().includes('name')
    ) || Object.keys(metrics[0])[0];

    const accuracyKey = Object.keys(metrics[0]).find(key =>
      key.toLowerCase().includes('accuracy') || key.toLowerCase().includes('acc')
    ) || Object.keys(metrics[0]).find(key => typeof metrics[0][key] === 'number');

    // Find best model based on accuracy
    const bestModel = [...metrics].sort((a, b) => {
      const aVal = typeof a[accuracyKey] === 'number' ? a[accuracyKey] : 0;
      const bVal = typeof b[accuracyKey] === 'number' ? b[accuracyKey] : 0;
      return bVal - aVal;
    })[0];

    // Calculate average accuracy
    const avgAccuracy = metrics.reduce((sum, m) => {
      const val = typeof m[accuracyKey] === 'number' ? m[accuracyKey] : 0;
      return sum + val;
    }, 0) / metrics.length;

    return {
      bestModel: bestModel[modelNameKey] || "Model " + metrics.indexOf(bestModel),
      bestAccuracy: typeof bestModel[accuracyKey] === 'number' ? bestModel[accuracyKey] : 0,
      avgAccuracy: avgAccuracy,
      totalModels: metrics.length
    };
  };

  const summary = calculateSummary();

  if (loading) {
    return (
      <div className="flex items-center justify-center h-screen bg-gray-50">
        <div className="text-center">
          <div className="inline-block animate-spin rounded-full border-4 border-gray-300 border-t-blue-600 h-12 w-12 mb-4"></div>
          <h2 className="text-xl font-semibold text-gray-700">Loading classification results...</h2>
          <p className="text-gray-500 mt-2">Please wait while we process the data</p>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center h-screen bg-gray-50">
        <div className="max-w-md w-full bg-white p-8 rounded-lg shadow-lg">
          <div className="flex items-center justify-center h-12 w-12 rounded-full bg-red-100 mb-4 mx-auto">
            <AlertTriangle className="h-6 w-6 text-red-600" />
          </div>
          <h2 className="text-center text-2xl font-bold text-gray-800 mb-2">Error Loading Data</h2>
          <p className="text-center text-gray-600 mb-6">{error}</p>
          <button
            className="w-full py-2 px-4 border border-transparent rounded-md shadow-sm text-white bg-blue-600 hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-blue-500"
            onClick={() => window.location.reload()}
          >
            Retry
          </button>
        </div>
      </div>
    );
  }

  const visibleMetricColumns = Object.keys(metrics[0] || {}).filter(col => visibleColumns[col]);
  const visiblePredictionColumns = Object.keys(predictions[0] || {});

  return (
    <div className="min-h-screen bg-gray-50 text-gray-800 font-mono">
      {/* Header */}
      <header className="bg-gray-800 text-white py-4 px-6 shadow-md">
        <div className="flex items-center justify-between">
          <div className="flex items-center space-x-2">
            <Database className="h-6 w-6" />
            <h1 className="text-xl font-bold tracking-tight">Clinical Classification Analysis</h1>
          </div>
          <div className="text-xs text-gray-400">
            Data fetched from: localhost:8000/clinical/classification
          </div>
        </div>
      </header>

      {/* Summary Stats */}
      <div className="container mx-auto p-6">
        <div className="grid grid-cols-1 md:grid-cols-4 gap-4 mb-6">
          <div className="bg-white rounded-lg shadow-md p-4 border-l-4 border-blue-500">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-xs text-gray-500 uppercase tracking-wider">Best Model</p>
                <p className="text-xl font-bold">{summary.bestModel}</p>
              </div>
              <CheckCircle className="h-8 w-8 text-blue-500" />
            </div>
          </div>

          <div className="bg-white rounded-lg shadow-md p-4 border-l-4 border-green-500">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-xs text-gray-500 uppercase tracking-wider">Best Accuracy</p>
                <p className="text-xl font-bold">{(summary.bestAccuracy * 100).toFixed(2)}%</p>
              </div>
              <ArrowUpRight className="h-8 w-8 text-green-500" />
            </div>
          </div>

          <div className="bg-white rounded-lg shadow-md p-4 border-l-4 border-yellow-500">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-xs text-gray-500 uppercase tracking-wider">Avg Accuracy</p>
                <p className="text-xl font-bold">{(summary.avgAccuracy * 100).toFixed(2)}%</p>
              </div>
              <BarChart2 className="h-8 w-8 text-yellow-500" />
            </div>
          </div>

          <div className="bg-white rounded-lg shadow-md p-4 border-l-4 border-purple-500">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-xs text-gray-500 uppercase tracking-wider">Models Tested</p>
                <p className="text-xl font-bold">{summary.totalModels}</p>
              </div>
              <FileText className="h-8 w-8 text-purple-500" />
            </div>
          </div>
        </div>

        {/* Tabs */}
        <div className="bg-white rounded-lg shadow-md mb-6">
          <div className="flex border-b border-gray-200">
            <button
              className={`px-6 py-3 font-medium text-sm ${activeTab === 'metrics' ? 'text-blue-600 border-b-2 border-blue-500' : 'text-gray-500 hover:text-gray-700'}`}
              onClick={() => setActiveTab('metrics')}
            >
              Performance Metrics
            </button>
            <button
              className={`px-6 py-3 font-medium text-sm ${activeTab === 'predictions' ? 'text-blue-600 border-b-2 border-blue-500' : 'text-gray-500 hover:text-gray-700'}`}
              onClick={() => setActiveTab('predictions')}
            >
              Holdout Predictions
            </button>
            <button
              className={`px-6 py-3 font-medium text-sm ${activeTab === 'visualization' ? 'text-blue-600 border-b-2 border-blue-500' : 'text-gray-500 hover:text-gray-700'}`}
              onClick={() => setActiveTab('visualization')}
            >
              Data Visualization
            </button>
          </div>
        </div>

        {/* Content based on active tab */}
        {activeTab === 'metrics' && (
          <div className="bg-white rounded-lg shadow-md p-6">
            <div className="flex justify-between items-center mb-4">
              <h2 className="text-lg font-bold">Model Performance Metrics</h2>
              <div className="flex space-x-2">
                <div className="relative">
                  <select
                    className="appearance-none bg-gray-50 border border-gray-300 text-gray-700 py-2 px-4 pr-8 rounded leading-tight focus:outline-none focus:bg-white focus:border-gray-500 text-sm"
                  >
                    <option>Sort by Accuracy</option>
                    <option>Sort by Precision</option>
                    <option>Sort by Recall</option>
                    <option>Sort by F1-Score</option>
                  </select>
                  <div className="pointer-events-none absolute inset-y-0 right-0 flex items-center px-2 text-gray-700">
                    <svg className="fill-current h-4 w-4" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 20 20">
                      <path d="M9.293 12.95l.707.707L15.657 8l-1.414-1.414L10 10.828 5.757 6.586 4.343 8z"/>
                    </svg>
                  </div>
                </div>
                <div className="relative">
                  <select
                    className="appearance-none bg-gray-50 border border-gray-300 text-gray-700 py-2 px-4 pr-8 rounded leading-tight focus:outline-none focus:bg-white focus:border-gray-500 text-sm"
                  >
                    <option>Filter Results</option>
                    <option>Accuracy > 0.8</option>
                    <option>Precision > 0.8</option>
                    <option>Recall > 0.8</option>
                  </select>
                  <div className="pointer-events-none absolute inset-y-0 right-0 flex items-center px-2 text-gray-700">
                    <svg className="fill-current h-4 w-4" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 20 20">
                      <path d="M9.293 12.95l.707.707L15.657 8l-1.414-1.414L10 10.828 5.757 6.586 4.343 8z"/>
                    </svg>
                  </div>
                </div>
              </div>
            </div>

            {/* Column selector */}
            <div className="flex flex-wrap gap-2 mb-4">
              {Object.keys(metrics[0] || {}).map(column => (
                <button
                  key={column}
                  onClick={() => toggleColumnVisibility(column)}
                  className={`px-2 py-1 rounded text-xs font-mono ${
                    visibleColumns[column] 
                      ? 'bg-blue-100 text-blue-800 border border-blue-200' 
                      : 'bg-gray-100 text-gray-500 border border-gray-200'
                  }`}
                >
                  {column}
                </button>
              ))}
            </div>

            {/* Metrics Table */}
            <div className="overflow-x-auto">
              <table className="min-w-full divide-y divide-gray-200">
                <thead className="bg-gray-50">
                  <tr>
                    {visibleMetricColumns.map(column => (
                      <th
                        key={column}
                        scope="col"
                        className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider"
                      >
                        {column}
                      </th>
                    ))}
                  </tr>
                </thead>
                <tbody className="bg-white divide-y divide-gray-200">
                  {metrics.map((row, rowIndex) => (
                    <tr key={rowIndex} className={rowIndex % 2 === 0 ? 'bg-gray-50' : 'bg-white'}>
                      {visibleMetricColumns.map(column => (
                        <td
                          key={column}
                          className="px-6 py-4 whitespace-nowrap text-sm font-mono"
                        >
                          {typeof row[column] === 'number'
                            ? row[column].toFixed(4)
                            : row[column]?.toString() || '-'}
                        </td>
                      ))}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}

        {activeTab === 'predictions' && (
          <div className="bg-white rounded-lg shadow-md p-6">
            <div className="flex justify-between items-center mb-4">
              <h2 className="text-lg font-bold">Holdout Predictions</h2>
              <div className="text-sm text-gray-500">
                Total predictions: {predictions.length}
              </div>
            </div>

            {/* Predictions Table */}
            <div className="overflow-x-auto" style={{ maxHeight: "60vh" }}>
              <table className="min-w-full divide-y divide-gray-200">
                <thead className="bg-gray-50 sticky top-0">
                  <tr>
                    {visiblePredictionColumns.map(column => (
                      <th
                        key={column}
                        scope="col"
                        className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider"
                      >
                        {column}
                      </th>
                    ))}
                  </tr>
                </thead>
                <tbody className="bg-white divide-y divide-gray-200">
                  {predictions.map((row, rowIndex) => (
                    <tr
                      key={rowIndex}
                      className={`${rowIndex % 2 === 0 ? 'bg-gray-50' : 'bg-white'} ${
                        Object.keys(row).some(k => k.toLowerCase() === 'actual') && 
                        Object.keys(row).some(k => k.toLowerCase() === 'predicted') &&
                        row[Object.keys(row).find(k => k.toLowerCase() === 'actual')] === 
                        row[Object.keys(row).find(k => k.toLowerCase() === 'predicted')] 
                          ? 'border-l-4 border-green-200' 
                          : Object.keys(row).some(k => k.toLowerCase() === 'actual') && 
                            Object.keys(row).some(k => k.toLowerCase() === 'predicted')
                              ? 'border-l-4 border-red-200'
                              : ''
                      }`}
                    >
                      {visiblePredictionColumns.map(column => (
                        <td
                          key={column}
                          className={`px-6 py-4 whitespace-nowrap text-sm font-mono ${
                            column.toLowerCase() === 'predicted' && 
                            Object.keys(row).some(k => k.toLowerCase() === 'actual') && 
                            row[column] !== row[Object.keys(row).find(k => k.toLowerCase() === 'actual')]
                              ? 'text-red-600 font-semibold' 
                              : column.toLowerCase() === 'predicted' && 
                                Object.keys(row).some(k => k.toLowerCase() === 'actual') && 
                                row[column] === row[Object.keys(row).find(k => k.toLowerCase() === 'actual')]
                                ? 'text-green-600 font-semibold'
                                : ''
                          }`}
                        >
                          {typeof row[column] === 'number'
                            ? row[column].toFixed(4)
                            : row[column]?.toString() || '-'}
                        </td>
                      ))}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}

        {activeTab === 'visualization' && (
          <div className="bg-white rounded-lg shadow-md p-6">
            <h2 className="text-lg font-bold mb-4">Performance Visualization</h2>

            <div className="h-96">
              <ResponsiveContainer width="100%" height="100%">
                <LineChart
                  data={chartData}
                  margin={{
                    top: 5,
                    right: 30,
                    left: 20,
                    bottom: 5,
                  }}
                >
                  <CartesianGrid strokeDasharray="3 3" />
                  <XAxis dataKey="name" />
                  <YAxis />
                  <Tooltip contentStyle={{ fontFamily: 'monospace' }} />
                  <Legend />
                  <Line type="monotone" dataKey="accuracy" stroke="#3B82F6" strokeWidth={2} activeDot={{ r: 8 }} />
                  <Line type="monotone" dataKey="precision" stroke="#10B981" strokeWidth={2} />
                  <Line type="monotone" dataKey="recall" stroke="#F59E0B" strokeWidth={2} />
                  <Line type="monotone" dataKey="f1" stroke="#8B5CF6" strokeWidth={2} />
                </LineChart>
              </ResponsiveContainer>
            </div>

            <div className="mt-6 bg-gray-50 border border-gray-200 rounded-md p-4">
              <h3 className="text-sm font-semibold mb-2">Analysis Notes</h3>
              <p className="text-sm text-gray-600">
                The chart above shows the performance metrics across different models.
                Higher values indicate better performance. The best performing model
                is <span className="font-semibold">{summary.bestModel}</span> with an
                accuracy of <span className="font-semibold">{(summary.bestAccuracy * 100).toFixed(2)}%</span>.
              </p>
            </div>
          </div>
        )}
      </div>

      {/* Footer */}
      <footer className="bg-gray-800 text-gray-400 py-4 px-6 text-center text-xs">
        <p>Classification Analytics Dashboard • {new Date().toLocaleDateString()}</p>
      </footer>
    </div>
  );
};

export default ClassificationDemo;