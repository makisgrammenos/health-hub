'use client';

import { useState, useEffect } from "react";
import { Database, AlertTriangle, BarChart2, Server, Activity, FileText, TrendingUp, Grid, Target } from "lucide-react";
import { ResponsiveContainer, LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend } from 'recharts';

const ClusteringVisualization = () => {
  const [clusterAssignments, setClusterAssignments] = useState([]);
  const [metrics, setMetrics] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [activeTab, setActiveTab] = useState('overview');
  const [selectedPlot, setSelectedPlot] = useState('pca_2d');


  const [modalOpen, setModalOpen] = useState(false);
  const [modalSrc, setModalSrc] = useState('');

  useEffect(() => {
    const fetchClusteringData = async () => {
      try {
        const response = await fetch("http://localhost:7000/clinical/clustering");
        if (!response.ok) throw new Error("Failed to fetch clustering data");
        const data = await response.json();
        setClusterAssignments(data.cluster_assignments);
        setMetrics(data.metrics);
      } catch (error) {
        console.error("Error fetching clustering data:", error);
        setError(error.message);
      } finally {
        setLoading(false);
      }
    };

    fetchClusteringData();
  }, []);

  // Calculate summary statistics
  const calculateSummary = () => {
    if (!clusterAssignments.length || !metrics.length) {
      return {
        totalClusters: 0,
        totalSamples: 0,
        silhouetteScore: 0,
        algorithm: "N/A"
      };
    }

    // Find number of unique clusters
    const clusterColumn = Object.keys(clusterAssignments[0]).find(key =>
      key.toLowerCase().includes('cluster') || key.toLowerCase().includes('label') || key.toLowerCase().includes('group')
    );

    const uniqueClusters = new Set();
    if (clusterColumn) {
      clusterAssignments.forEach(item => uniqueClusters.add(item[clusterColumn]));
    }

    // Find silhouette score if available
    const silhouetteKey = Object.keys(metrics[0] || {}).find(key =>
      key.toLowerCase().includes('silhouette') || key.toLowerCase().includes('score')
    );

    const silhouetteScore = silhouetteKey && metrics.length > 0 ?
      metrics.find(m => typeof m[silhouetteKey] === 'number')?.[silhouetteKey] || 0 : 0;

    // Find algorithm if available
    const algorithmKey = Object.keys(metrics[0] || {}).find(key =>
      key.toLowerCase().includes('algorithm') || key.toLowerCase().includes('method')
    );

    const algorithm = algorithmKey && metrics.length > 0 ?
      metrics.find(m => m[algorithmKey])?.[algorithmKey] || "KMeans" : "KMeans";

    return {
      totalClusters: uniqueClusters.size || 0,
      totalSamples: clusterAssignments.length,
      silhouetteScore: silhouetteScore,
      algorithm: algorithm
    };
  };

  const summary = calculateSummary();

  // Prepare distribution data for chart
  const prepareDistributionData = () => {
    if (!clusterAssignments.length) return [];

    const clusterColumn = Object.keys(clusterAssignments[0]).find(key =>
      key.toLowerCase().includes('cluster') || key.toLowerCase().includes('label') || key.toLowerCase().includes('group')
    );

    if (!clusterColumn) return [];

    // Count items per cluster
    const counts = {};
    clusterAssignments.forEach(item => {
      const cluster = item[clusterColumn];
      counts[cluster] = (counts[cluster] || 0) + 1;
    });

    // Convert to chart data
    return Object.entries(counts).map(([cluster, count]) => ({
      name: `Cluster ${cluster}`,
      count: count
    }));
  };
  const openModal = (src: string) => {
    setModalSrc(src);
    setModalOpen(true);
  };
  const closeModal = () => setModalOpen(false);

  const distributionData = prepareDistributionData();

  if (loading) {
    return (
      <div className="flex items-center justify-center h-screen bg-gray-50">
        <div className="text-center">
          <div className="inline-block animate-spin rounded-full border-4 border-gray-300 border-t-blue-600 h-12 w-12 mb-4"></div>
          <h2 className="text-xl font-semibold text-gray-700">Loading clustering results...</h2>
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

  const plotImages = [
    { id: 'pca_2d', name: 'PCA 2D Plot', src: '/clustering/pca_2d.png', icon: Grid },
    { id: 'pca_2d_tag', name: 'PCA with Labels', src: '/clustering/pca_2d_tag.png', icon: Target },
    { id: 'tsne_3d', name: '3D T-SNE Plot', src: '/clustering/tsne_3d.png', icon: TrendingUp },
    { id: 'silhouette', name: 'Silhouette Plot', src: '/clustering/silhouette.png', icon: Activity },
    { id: 'kmeans_intercluster', name: 'Distance Plot', src: '/clustering/kmeans_intercluster.png', icon: Server },
    { id: 'cluster_dist', name: 'Distribution Plot', src: '/clustering/cluster_dist.png', icon: BarChart2 }
  ];

  return (
    <div className="min-h-screen bg-gray-50 text-gray-800 font-mono">
      {/* Header */}
      <header className="bg-gray-800 text-white py-4 px-6 shadow-md">
        <div className="flex items-center justify-between">
          <div className="flex items-center space-x-2">
            <Database className="h-6 w-6" />
            <h1 className="text-xl font-bold tracking-tight">Clustering Analysis Dashboard</h1>
          </div>
          <div className="text-xs text-gray-400">
            Data fetched from: localhost:7000/clinical/clustering
          </div>
        </div>
      </header>

      {/* Summary Stats */}
      <div className="container mx-auto p-6">
        <div className="grid grid-cols-1 md:grid-cols-4 gap-4 mb-6">
          <div className="bg-white rounded-lg shadow-md p-4 border-l-4 border-blue-500">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-xs text-gray-500 uppercase tracking-wider">Algorithm</p>
                <p className="text-xl font-bold">{summary.algorithm}</p>
              </div>
              <FileText className="h-8 w-8 text-blue-500" />
            </div>
          </div>

          <div className="bg-white rounded-lg shadow-md p-4 border-l-4 border-green-500">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-xs text-gray-500 uppercase tracking-wider">Total Clusters</p>
                <p className="text-xl font-bold">{summary.totalClusters}</p>
              </div>
              <Target className="h-8 w-8 text-green-500" />
            </div>
          </div>

          <div className="bg-white rounded-lg shadow-md p-4 border-l-4 border-yellow-500">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-xs text-gray-500 uppercase tracking-wider">Total Samples</p>
                <p className="text-xl font-bold">{summary.totalSamples}</p>
              </div>
              <Database className="h-8 w-8 text-yellow-500" />
            </div>
          </div>

          <div className="bg-white rounded-lg shadow-md p-4 border-l-4 border-purple-500">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-xs text-gray-500 uppercase tracking-wider">Silhouette Score</p>
                <p className="text-xl font-bold">{summary.silhouetteScore.toFixed(4)}</p>
              </div>
              <Activity className="h-8 w-8 text-purple-500" />
            </div>
          </div>
        </div>

        {/* Tabs */}
        <div className="bg-white rounded-lg shadow-md mb-6">
          <div className="flex border-b border-gray-200">
            <button
              className={`px-6 py-3 font-medium text-sm ${activeTab === 'overview' ? 'text-blue-600 border-b-2 border-blue-500' : 'text-gray-500 hover:text-gray-700'}`}
              onClick={() => setActiveTab('overview')}
            >
              Overview
            </button>
            <button
              className={`px-6 py-3 font-medium text-sm ${activeTab === 'visualizations' ? 'text-blue-600 border-b-2 border-blue-500' : 'text-gray-500 hover:text-gray-700'}`}
              onClick={() => setActiveTab('visualizations')}
            >
              Visualizations
            </button>
            <button
              className={`px-6 py-3 font-medium text-sm ${activeTab === 'data' ? 'text-blue-600 border-b-2 border-blue-500' : 'text-gray-500 hover:text-gray-700'}`}
              onClick={() => setActiveTab('data')}
            >
              Raw Data
            </button>
            <button
              className={`px-6 py-3 font-medium text-sm ${activeTab === 'metrics' ? 'text-blue-600 border-b-2 border-blue-500' : 'text-gray-500 hover:text-gray-700'}`}
              onClick={() => setActiveTab('metrics')}
            >
              Metrics
            </button>
          </div>
        </div>

        {/* Content based on active tab */}
        {activeTab === 'overview' && (
          <div className="space-y-6">
            <div className="bg-white rounded-lg shadow-md p-6">
              <h2 className="text-lg font-bold mb-4">Cluster Distribution</h2>
              <div className="h-64">
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart
                    data={distributionData}
                    margin={{ top: 5, right: 30, left: 20, bottom: 20 }}
                  >
                    <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" />
                    <XAxis dataKey="name" tick={{ fontSize: 12 }} />
                    <YAxis tick={{ fontSize: 12 }} />
                    <Tooltip contentStyle={{ fontFamily: 'monospace', fontSize: 12 }} />
                    <Legend />
                    <Line
                      type="monotone"
                      dataKey="count"
                      name="Sample Count"
                      stroke="#3B82F6"
                      strokeWidth={3}
                      dot={{ r: 4 }}
                      activeDot={{ r: 6 }}
                    />
                  </LineChart>
                </ResponsiveContainer>
              </div>
            </div>

            <div className="bg-white rounded-lg shadow-md p-6">
              <h2 className="text-lg font-bold mb-4">PCA Visualization</h2>
              <img
                src="/clustering/pca_2d.png"
                alt="PCA Cluster Plot"
                className="w-full max-h-96 object-contain rounded border border-gray-200"
                onClick={()=>openModal(<img
                src="/clustering/pca_2d.png"
                alt="PCA Cluster Plot"
                className="w-full max-h-96 object-contain rounded border border-gray-200"

              />)}
              />
              <div className="mt-4 p-4 bg-gray-50 border border-gray-200 rounded-md">
                <h3 className="text-sm font-semibold mb-2">Analysis Notes</h3>
                <p className="text-sm text-gray-700">
                  Principal Component Analysis (PCA) reduces data dimensionality while preserving variance.
                  This 2D projection shows {summary.totalClusters} distinct clusters identified by the {summary.algorithm} algorithm.
                  The spatial separation between clusters indicates the quality of the clustering solution.
                </p>
              </div>
            </div>

            <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
              <div className="bg-white rounded-lg shadow-md p-6">
                <h2 className="text-lg font-bold mb-4">Silhouette Analysis</h2>
                <img
                  src="/clustering/silhouette.png"
                  alt="Silhouette Plot"
                  className="w-full max-h-64 object-contain rounded border border-gray-200"
                />
                <div className="mt-4 p-3 bg-gray-50 border border-gray-200 rounded-md">
                  <p className="text-xs text-gray-700">
                    Silhouette score: <span className="font-bold">{summary.silhouetteScore.toFixed(4)}</span>.
                    Values closer to 1 indicate better-defined clusters.
                  </p>
                </div>
              </div>

              <div className="bg-white rounded-lg shadow-md p-6">
                <h2 className="text-lg font-bold mb-4">Cluster Distribution</h2>
                <img
                  src="/clustering/cluster_dist.png"
                  alt="Distribution Plot"
                  className="w-full max-h-64 object-contain rounded border border-gray-200"
                />
                <div className="mt-4 p-3 bg-gray-50 border border-gray-200 rounded-md">
                  <p className="text-xs text-gray-700">
                    Distribution of {summary.totalSamples} samples across {summary.totalClusters} clusters,
                    showing the relative cluster sizes and balance.
                  </p>
                </div>
              </div>
            </div>
          </div>
        )}

        {activeTab === 'visualizations' && (
          <div className="bg-white rounded-lg shadow-md p-6">
            <div className="flex justify-between items-center mb-6">
              <h2 className="text-lg font-bold">Cluster Visualizations</h2>
              <div className="text-sm text-gray-500">
                {summary.totalClusters} clusters · {summary.totalSamples} samples
              </div>
            </div>

            <div className="grid grid-cols-2 md:grid-cols-3 gap-4 mb-6">
              {plotImages.map((plot) => (
                <button
                  key={plot.id}
                  onClick={() => setSelectedPlot(plot.id)}
                  className={`flex items-center p-3 rounded-md border ${
                    selectedPlot === plot.id 
                      ? 'border-blue-500 bg-blue-50 text-blue-700' 
                      : 'border-gray-200 hover:bg-gray-50'
                  }`}
                >
                  <plot.icon className="h-5 w-5 mr-2" />
                  <span className="text-sm font-medium">{plot.name}</span>
                </button>
              ))}
            </div>

            <div className="bg-gray-50 p-1 rounded-lg border border-gray-200">
              <img
                src={plotImages.find(p => p.id === selectedPlot)?.src || '/clustering/pca_2d.png'}
                alt={plotImages.find(p => p.id === selectedPlot)?.name || 'Cluster Visualization'}
                className="w-full h-auto max-h-[28rem] object-contain rounded mx-auto"
              />
            </div>

            <div className="mt-6 p-4 bg-gray-50 border border-gray-200 rounded-md">
              <h3 className="text-sm font-semibold mb-2">Visualization Description</h3>

              {selectedPlot === 'pca_2d' && (
                <p className="text-sm text-gray-700">
                  This 2D PCA plot shows the data points projected onto the first two principal components.
                  Different colors represent different clusters identified by the {summary.algorithm} algorithm.
                  The spatial separation indicates how well-defined the clusters are in the reduced feature space.
                </p>
              )}

              {selectedPlot === 'pca_2d_tag' && (
                <p className="text-sm text-gray-700">
                  PCA projection with labeled data points. This visualization shows how the original class labels
                  relate to the discovered clusters, allowing for evaluation of clustering accuracy when ground truth is available.
                </p>
              )}

              {selectedPlot === 'tsne_3d' && (
                <p className="text-sm text-gray-700">
                  t-SNE (t-Distributed Stochastic Neighbor Embedding) creates a 3D representation that preserves
                  local similarities in the high-dimensional space. This non-linear technique often reveals more
                  complex patterns than PCA.
                </p>
              )}

              {selectedPlot === 'silhouette' && (
                <p className="text-sm text-gray-700">
                  The silhouette plot displays how well each point lies within its cluster.
                  Values near +1 indicate points well-matched to their clusters and far from neighboring clusters.
                  The average silhouette width of {summary.silhouetteScore.toFixed(4)} indicates the overall clustering quality.
                </p>
              )}

              {selectedPlot === 'kmeans_intercluster' && (
                <p className="text-sm text-gray-700">
                  The distance plot visualizes the distances between cluster centroids.
                  Larger distances between centroids generally indicate better cluster separation.
                  This helps evaluate the distinctiveness of the identified clusters.
                </p>
              )}

              {selectedPlot === 'cluster_dist' && (
                <p className="text-sm text-gray-700">
                  This distribution plot shows how samples are distributed across the {summary.totalClusters} clusters.
                  Ideally, clusters should be relatively balanced, though some variation is expected due to
                  the natural structure of the data.
                </p>
              )}
            </div>
          </div>
        )}

        {activeTab === 'data' && (
          <div className="bg-white rounded-lg shadow-md p-6">
            <div className="flex justify-between items-center mb-4">
              <h2 className="text-lg font-bold">Cluster Assignments</h2>
              <div className="text-sm text-gray-500">
                {clusterAssignments.length} records
              </div>
            </div>

            <div className="overflow-x-auto" style={{ maxHeight: "60vh" }}>
              <table className="min-w-full divide-y divide-gray-200">
                <thead className="bg-gray-50 sticky top-0">
                  <tr>
                    {Object.keys(clusterAssignments[0] || {}).map((col) => (
                      <th
                        key={col}
                        scope="col"
                        className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider"
                      >
                        {col}
                      </th>
                    ))}
                  </tr>
                </thead>
                <tbody className="bg-white divide-y divide-gray-200">
                  {clusterAssignments.map((item, index) => {
                    const clusterKey = Object.keys(item).find(key =>
                      key.toLowerCase().includes('cluster') || key.toLowerCase().includes('label')
                    );

                    const clusterValue = clusterKey ? item[clusterKey] : null;

                    return (
                      <tr
                        key={index}
                        className={`${index % 2 === 0 ? 'bg-gray-50' : 'bg-white'} hover:bg-gray-100`}
                      >
                        {Object.keys(item).map((key) => (
                          <td
                            key={key}
                            className={`px-6 py-4 whitespace-nowrap text-sm font-mono ${
                              key === clusterKey ? 'font-bold text-blue-600' : ''
                            }`}
                          >
                            {typeof item[key] === 'number'
                              ? item[key].toFixed(4)
                              : item[key]?.toString() || '-'}
                          </td>
                        ))}
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>
          </div>
        )}

        {activeTab === 'metrics' && (
          <div className="bg-white rounded-lg shadow-md p-6">
            <div className="flex justify-between items-center mb-4">
              <h2 className="text-lg font-bold">Clustering Metrics</h2>
              <div className="text-sm text-gray-500">
                Evaluation results
              </div>
            </div>

            <div className="overflow-x-auto">
              <table className="min-w-full divide-y divide-gray-200">
                <thead className="bg-gray-50">
                  <tr>
                    {Object.keys(metrics[0] || {}).map((col) => (
                      <th
                        key={col}
                        scope="col"
                        className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider"
                      >
                        {col}
                      </th>
                    ))}
                  </tr>
                </thead>
                <tbody className="bg-white divide-y divide-gray-200">
                  {metrics.map((item, index) => (
                    <tr key={index} className={index % 2 === 0 ? 'bg-gray-50' : 'bg-white'}>
                      {Object.keys(item).map((key) => (
                        <td
                          key={key}
                          className="px-6 py-4 whitespace-nowrap text-sm font-mono"
                        >
                          {typeof item[key] === 'number'
                            ? item[key].toFixed(4)
                            : item[key]?.toString() || '-'}
                        </td>
                      ))}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>

            <div className="mt-6 p-4 bg-gray-50 border border-gray-200 rounded-md">
              <h3 className="text-sm font-semibold mb-2">Metrics Overview</h3>
              <p className="text-sm text-gray-700">
                The clustering evaluation metrics provide quantitative assessment of the clustering quality.
                Key metrics include:
              </p>
              <ul className="list-disc pl-5 mt-2 text-sm text-gray-700 space-y-1">
                <li>
                  <span className="font-medium">Silhouette Score:</span> Measures how similar objects are to their own cluster compared to other clusters (range: -1 to 1)
                </li>
                <li>
                  <span className="font-medium">Inertia:</span> Sum of squared distances to the nearest centroid (lower is better)
                </li>
                <li>
                  <span className="font-medium">Calinski-Harabasz Index:</span> Ratio of between-cluster to within-cluster dispersion (higher is better)
                </li>
                <li>
                  <span className="font-medium">Davies-Bouldin Index:</span> Average similarity between clusters (lower is better)
                </li>
              </ul>
            </div>
          </div>
        )}
      </div>

      {/* Footer */}
      <footer className="bg-gray-800 text-gray-400 py-4 px-6 text-center text-xs">
        <p>Clustering Analysis Dashboard • {new Date().toLocaleDateString()}</p>
      </footer>
    </div>
  );
};

export default ClusteringVisualization;