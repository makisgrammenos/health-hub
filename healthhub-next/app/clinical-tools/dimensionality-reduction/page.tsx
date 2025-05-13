'use client';

import { useState, useEffect } from "react";
import {
  Table,
  TableHeader,
  TableColumn,
  TableBody,
  TableRow,
  TableCell,
  Card,
  CardHeader,
  CardBody,
  Image,
  Tabs,
  Tab,
  Dropdown,
  DropdownTrigger,
  DropdownMenu,
  DropdownItem,
  Button,
  Progress,
  Tooltip,
  Chip
} from "@nextui-org/react";

const DimReductionVisualization = () => {
  const [pca2D, setPca2D] = useState([]);
  const [umap2D, setUmap2D] = useState([]);
  const [tsne2D, setTsne2D] = useState([]);
  const [pca3D, setPca3D] = useState([]);
  const [umap3D, setUmap3D] = useState([]);
  const [tsne3D, setTsne3D] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [activeTab, setActiveTab] = useState("2d");
  const [selectedAlgorithm, setSelectedAlgorithm] = useState("pca");
  const [loadingProgress, setLoadingProgress] = useState({
    pca2D: 0,
    umap2D: 0,
    tsne2D: 0,
    pca3D: 0,
    umap3D: 0,
    tsne3D: 0,
  });

  useEffect(() => {
    const fetchData = async () => {
      try {
        const fetchDataset = async (url, dataType) => {
          // Simulate loading progress for better UX
          const simulateProgress = () => {
            let progress = 0;
            const interval = setInterval(() => {
              progress += Math.floor(Math.random() * 15) + 5;
              setLoadingProgress(prev => ({
                ...prev,
                [dataType]: progress > 100 ? 100 : progress
              }));
              if (progress >= 100) clearInterval(interval);
            }, 300);
            return () => clearInterval(interval);
          };

          const cleanup = simulateProgress();
          const response = await fetch(url);
          if (!response.ok) throw new Error(`Failed to fetch data from ${url}`);
          cleanup();
          setLoadingProgress(prev => ({ ...prev, [dataType]: 100 }));
          return response.json();
        };

        // Fetch all datasets
        setPca2D(await fetchDataset("http://localhost:8000/clinical/dim_reduction/pca_2d", "pca2D"));
        setUmap2D(await fetchDataset("http://localhost:8000/clinical/dim_reduction/umap_2d", "umap2D"));
        setTsne2D(await fetchDataset("http://localhost:8000/clinical/dim_reduction/tsne_2d", "tsne2D"));
        setPca3D(await fetchDataset("http://localhost:8000/clinical/dim_reduction/pca_3d", "pca3D"));
        setUmap3D(await fetchDataset("http://localhost:8000/clinical/dim_reduction/umap_3d", "umap3D"));
        setTsne3D(await fetchDataset("http://localhost:8000/clinical/dim_reduction/tsne_3d", "tsne3D"));
      } catch (err) {
        console.error("Error fetching data:", err);
        setError(err.message);
      } finally {
        setLoading(false);
      }
    };

    fetchData();
  }, []);

  const algorithmInfo = {
    pca: {
      name: "Principal Component Analysis",
      description: "PCA reduces dimensionality by projecting data onto orthogonal axes that maximize variance. It preserves global structure but may lose local relationships.",
      citation: "Pearson, K. (1901). On lines and planes of closest fit to systems of points in space. Philosophical Magazine, 2(11), 559-572.",
      complexity: "O(min(n²d, nd²))",
      strengths: ["Linear technique", "Handles noise well", "Computationally efficient", "Provides variance explained"],
      weaknesses: ["Cannot capture non-linear relationships", "May lose important local structure", "Assumes orthogonality"]
    },
    umap: {
      name: "Uniform Manifold Approximation and Projection",
      description: "UMAP constructs a high-dimensional graph representation and optimizes a low-dimensional layout to preserve topological structure.",
      citation: "McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. ArXiv e-prints.",
      complexity: "O(n log n)",
      strengths: ["Preserves local and global structure", "Computationally efficient", "Works with large datasets", "Theoretical foundation"],
      weaknesses: ["Stochastic results", "Sensitive to hyperparameters", "Less interpretable than PCA"]
    },
    tsne: {
      name: "t-Distributed Stochastic Neighbor Embedding",
      description: "t-SNE converts similarities between data points to joint probabilities and minimizes the KL divergence between probability distributions.",
      citation: "van der Maaten, L., & Hinton, G. (2008). Visualizing Data using t-SNE. Journal of Machine Learning Research, 9, 2579-2605.",
      complexity: "O(n²)",
      strengths: ["Excellent at preserving local structure", "Handles non-linear relationships", "Reveals clusters effectively"],
      weaknesses: ["Computationally intensive", "Does not preserve global structure well", "Cannot meaningfully project new data"]
    }
  };

  const renderDataStatus = () => {
    const datasets = [
      { name: "PCA 2D", status: loadingProgress.pca2D },
      { name: "UMAP 2D", status: loadingProgress.umap2D },
      { name: "t-SNE 2D", status: loadingProgress.tsne2D },
      { name: "PCA 3D", status: loadingProgress.pca3D },
      { name: "UMAP 3D", status: loadingProgress.umap3D },
      { name: "t-SNE 3D", status: loadingProgress.tsne3D }
    ];

    return (
      <div className="bg-gray-100 p-4 rounded-lg mb-6">
        <h3 className="text-lg font-mono mb-3">Dataset Status</h3>
        <div className="grid grid-cols-2 md:grid-cols-3 gap-4">
          {datasets.map((dataset) => (
            <div key={dataset.name} className="mb-2">
              <div className="flex justify-between mb-1">
                <span className="text-sm font-mono">{dataset.name}</span>
                <span className="text-xs font-mono">{dataset.status}%</span>
              </div>
              <Progress
                color={dataset.status === 100 ? "success" : "primary"}
                value={dataset.status}
                size="sm"
              />
            </div>
          ))}
        </div>
      </div>
    );
  };

  const renderTable = (data, label) => (
    <div className="mt-6 rounded-lg border border-gray-200 overflow-hidden">
      <div className="bg-gray-50 border-b border-gray-200 py-2 px-4">
        <div className="flex justify-between items-center">
          <h3 className="font-mono text-sm">Numerical Data ({data.length} records)</h3>
          <Chip size="sm" variant="flat" color="primary" className="font-mono">
            {label}
          </Chip>
        </div>
      </div>
      <div className="max-h-64 overflow-y-auto">
        <Table
          aria-label={`${label} Data Table`}
          selectionMode="none"
          classNames={{
            base: "font-mono text-xs",
          }}
        >
          <TableHeader>
            {Object.keys(data[0] || {}).map((col) => (
              <TableColumn key={col}>{col}</TableColumn>
            ))}
          </TableHeader>
          <TableBody>
            {data.map((item, index) => (
              <TableRow key={index}>
                {Object.keys(item).map((key) => (
                  <TableCell key={key}>{item[key]}</TableCell>
                ))}
              </TableRow>
            ))}
          </TableBody>
        </Table>
      </div>
    </div>
  );

  const renderAlgorithmInfo = (algorithm) => {
    const info = algorithmInfo[algorithm];

    return (
      <div className="bg-gray-50 rounded-lg p-4 my-4 border border-gray-200">
        <h3 className="font-mono text-lg font-medium mb-2">{info.name}</h3>
        <p className="text-sm mb-3">{info.description}</p>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mt-4">
          <div>
            <h4 className="font-mono text-xs uppercase tracking-wider text-gray-500 mb-1">Computational Complexity</h4>
            <p className="font-mono bg-gray-100 p-2 rounded text-sm">{info.complexity}</p>
          </div>
          <div>
            <h4 className="font-mono text-xs uppercase tracking-wider text-gray-500 mb-1">Citation</h4>
            <p className="text-xs italic">{info.citation}</p>
          </div>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mt-4">
          <div>
            <h4 className="font-mono text-xs uppercase tracking-wider text-gray-500 mb-1">Strengths</h4>
            <ul className="list-disc list-inside text-xs">
              {info.strengths.map((strength, i) => (
                <li key={i}>{strength}</li>
              ))}
            </ul>
          </div>
          <div>
            <h4 className="font-mono text-xs uppercase tracking-wider text-gray-500 mb-1">Limitations</h4>
            <ul className="list-disc list-inside text-xs">
              {info.weaknesses.map((weakness, i) => (
                <li key={i}>{weakness}</li>
              ))}
            </ul>
          </div>
        </div>
      </div>
    );
  };

  const renderVisualization = () => {
    let data;
    let imagePath;
    let label;

    if (activeTab === "2d") {
      if (selectedAlgorithm === "pca") {
        data = pca2D;
        imagePath = "/dim/pca_2d.png";
        label = "PCA 2D";
      } else if (selectedAlgorithm === "umap") {
        data = umap2D;
        imagePath = "/dim/umap2d.png";
        label = "UMAP 2D";
      } else {
        data = tsne2D;
        imagePath = "/dim/tsne_2d.png";
        label = "t-SNE 2D";
      }
    } else {
      if (selectedAlgorithm === "pca") {
        data = pca3D;
        imagePath = "/dim/pca_3d.png";
        label = "PCA 3D";
      } else if (selectedAlgorithm === "umap") {
        data = umap3D;
        imagePath = "/dim/umap3d.png";
        label = "UMAP 3D";
      } else {
        data = tsne3D;
        imagePath = "/dim/tsne_3d.png";
        label = "t-SNE 3D";
      }
    }

    return (
      <div>
        {renderAlgorithmInfo(selectedAlgorithm)}
        <Card className="border border-gray-200">
          <CardHeader className="bg-gray-50 border-b border-gray-200">
            <div className="flex w-full justify-between items-center">
              <h2 className="font-mono text-lg">{label} Projection</h2>
              <div className="flex space-x-2">
                <Tooltip content="Download data as CSV">
                  <Button size="sm" isIconOnly variant="flat">
                    <DownloadIcon />
                  </Button>
                </Tooltip>
                <Tooltip content="Download image">
                  <Button size="sm" isIconOnly variant="flat">
                    <ImageIcon />
                  </Button>
                </Tooltip>
              </div>
            </div>
          </CardHeader>
          <CardBody>
            <div className="bg-black rounded-lg p-2 mb-4">
              <div className="relative pt-[100.00%]">
                <div className="absolute inset-0 flex items-center justify-center">
                  <Image
                    src={imagePath}
                    alt={`${label} Visualization`}
                    className="object-contain w-full h-full"
                  />
                </div>
                <div className="absolute top-2 right-2">
                  <Chip size="sm" variant="solid" color="primary" className="font-mono bg-opacity-70">
                    {label}
                  </Chip>
                </div>
                <div className="absolute bottom-2 left-2 bg-black bg-opacity-40 text-white p-1 rounded text-xs font-mono">
                  n = {data?.length || 0} points
                </div>
              </div>
            </div>
            {data?.length > 0 && renderTable(data, label)}
          </CardBody>
        </Card>
      </div>
    );
  };

  if (loading) {
    return (
      <div className="min-h-screen flex flex-col items-center justify-center p-4 bg-gray-50">
        <div className="w-full max-w-3xl">
          <h1 className="text-2xl font-mono text-center mb-8">
            Loading Dimensional Reduction Data...
          </h1>
          {renderDataStatus()}
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="min-h-screen flex items-center justify-center p-4 bg-gray-50">
        <Card className="w-full max-w-3xl border-red-300">
          <CardHeader className="bg-red-50 text-red-700">
            <h2 className="font-mono">Error Loading Data</h2>
          </CardHeader>
          <CardBody>
            <p className="font-mono text-sm">{error}</p>
            <Button color="primary" className="mt-4">
              Retry Loading Data
            </Button>
          </CardBody>
        </Card>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gray-50 font-sans">
      <div className="bg-gradient-to-r from-blue-900 to-purple-900 text-white py-6 px-6 mb-8">
        <h1 className="text-center font-mono text-3xl mb-2 font-bold tracking-tight">
          Multidimensional Data Analysis
        </h1>
        <p className="text-center text-gray-200 font-light">
          Visualizing high-dimensional data through advanced dimensionality reduction techniques
        </p>
      </div>

      <div className="container mx-auto px-4 pb-16 max-w-6xl">
        <div className="flex flex-col md:flex-row gap-6 mb-8">
          <Card className="md:w-1/3 border border-gray-200">
            <CardHeader className="bg-gray-50 border-b border-gray-200">
              <h2 className="font-mono">Control Panel</h2>
            </CardHeader>
            <CardBody>
              <div className="mb-6">
                <h3 className="text-sm font-mono mb-2">Projection Dimensions</h3>
                <Tabs
                  aria-label="Dimensions"
                  selectedKey={activeTab}
                  onSelectionChange={setActiveTab}
                  variant="bordered"
                  className="font-mono"
                >
                  <Tab key="2d" title="2D Projection" />
                  <Tab key="3d" title="3D Projection" />
                </Tabs>
              </div>

              <div>
                <h3 className="text-sm font-mono mb-2">Algorithm Selection</h3>
                <div className="grid grid-cols-1 gap-2">
                  <Button
                    variant={selectedAlgorithm === "pca" ? "solid" : "flat"}
                    color={selectedAlgorithm === "pca" ? "primary" : "default"}
                    className="font-mono justify-start"
                    onClick={() => setSelectedAlgorithm("pca")}
                  >
                    PCA
                  </Button>
                  <Button
                    variant={selectedAlgorithm === "umap" ? "solid" : "flat"}
                    color={selectedAlgorithm === "umap" ? "primary" : "default"}
                    className="font-mono justify-start"
                    onClick={() => setSelectedAlgorithm("umap")}
                  >
                    UMAP
                  </Button>
                  <Button
                    variant={selectedAlgorithm === "tsne" ? "solid" : "flat"}
                    color={selectedAlgorithm === "tsne" ? "primary" : "default"}
                    className="font-mono justify-start"
                    onClick={() => setSelectedAlgorithm("tsne")}
                  >
                    t-SNE
                  </Button>
                </div>
              </div>
            </CardBody>
          </Card>

          <div className="md:w-2/3">
            {renderVisualization()}
          </div>
        </div>

        <Card className="border border-gray-200 mb-8">
          <CardHeader className="bg-gray-50 border-b border-gray-200">
            <h2 className="font-mono">Methodology Overview</h2>
          </CardHeader>
          <CardBody>
            <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
              <div className="border-l-2 border-blue-500 pl-4">
                <h3 className="font-mono text-blue-700 mb-1">Data Pre-processing</h3>
                <p className="text-sm">The original high-dimensional dataset undergoes standardization (z-score normalization) followed by outlier detection and handling before dimensionality reduction techniques are applied.</p>
              </div>
              <div className="border-l-2 border-purple-500 pl-4">
                <h3 className="font-mono text-purple-700 mb-1">Dimensionality Reduction</h3>
                <p className="text-sm">Multiple techniques are employed to reduce the high-dimensional feature space to 2D and 3D representations while preserving the underlying structure and relationships.</p>
              </div>
              <div className="border-l-2 border-green-500 pl-4">
                <h3 className="font-mono text-green-700 mb-1">Validation &amp; Comparison</h3>
                <p className="text-sm">Results from different methods are quantitatively compared using metrics such as reconstruction error, trustworthiness, and continuity to assess performance.</p>
              </div>
            </div>
          </CardBody>
        </Card>

        <div className="text-center mt-12 text-sm text-gray-500">
          <p>Research Visualization Tool v1.2.0 | Data Processing Pipeline | Clinical Data Analysis Platform</p>
        </div>
      </div>
    </div>
  );
};

// Simple icon components
const DownloadIcon = () => (
  <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" viewBox="0 0 16 16">
    <path d="M.5 9.9a.5.5 0 0 1 .5.5v2.5a1 1 0 0 0 1 1h12a1 1 0 0 0 1-1v-2.5a.5.5 0 0 1 1 0v2.5a2 2 0 0 1-2 2H2a2 2 0 0 1-2-2v-2.5a.5.5 0 0 1 .5-.5z"/>
    <path d="M7.646 11.854a.5.5 0 0 0 .708 0l3-3a.5.5 0 0 0-.708-.708L8.5 10.293V1.5a.5.5 0 0 0-1 0v8.793L5.354 8.146a.5.5 0 1 0-.708.708l3 3z"/>
  </svg>
);

const ImageIcon = () => (
  <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" viewBox="0 0 16 16">
    <path d="M6.002 5.5a1.5 1.5 0 1 1-3 0 1.5 1.5 0 0 1 3 0z"/>
    <path d="M2.002 1a2 2 0 0 0-2 2v10a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V3a2 2 0 0 0-2-2h-12zm12 1a1 1 0 0 1 1 1v6.5l-3.777-1.947a.5.5 0 0 0-.577.093l-3.71 3.71-2.66-1.772a.5.5 0 0 0-.63.062L1.002 12V3a1 1 0 0 1 1-1h12z"/>
  </svg>
);

export default DimReductionVisualization;