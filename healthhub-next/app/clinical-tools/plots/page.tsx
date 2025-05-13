'use client';

import { useState, useEffect } from "react";
import {
  Table,
  TableHeader,
  TableColumn,
  TableBody,
  TableRow,
  TableCell,
} from "@nextui-org/react";
import { Card, CardBody, CardHeader, CardFooter } from "@nextui-org/card";
import { Divider } from "@nextui-org/divider";
import { Button } from "@nextui-org/button";
import { Tooltip } from "@nextui-org/tooltip";
import { Progress } from "@nextui-org/progress";
import { Chip } from "@nextui-org/chip";
import { Tabs, Tab } from "@nextui-org/tabs";

// Icons
import {
  BarChart2,
  Database,
  FileSpreadsheet,
  Grid,
  Layers,
  LineChart,
  Loader2,
  AlertTriangle,
  Download,
  Share2,
  Info,
  Maximize2,
  PieChart,
} from "lucide-react";

const DataVisualization = () => {
  const [categoryCounts, setCategoryCounts] = useState(null);
  const [standardizedData, setStandardizedData] = useState([]);
  const [meltedData, setMeltedData] = useState([]);
  const [correlationMatrix, setCorrelationMatrix] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [activeTab, setActiveTab] = useState("overview");
  const [fullscreenImage, setFullscreenImage] = useState(null);

  useEffect(() => {
    const fetchData = async () => {
      try {
        // Fetch category counts
       // const categoryResponse = await fetch("http://localhost:8000/clinical/heart/category_counts");
       //  if (!categoryResponse.ok) throw new Error("Failed to fetch category counts");
       //  const categoryJson = await categoryResponse.json();
       //  console.log("Category counts data:", categoryJson); // Debug
       //  setCategoryCounts(categoryJson.category_counts || categoryJson);
          const categoryResponse = await fetch("http://localhost:8000/clinical/heart/category_counts");
          if (!categoryResponse.ok) throw new Error("Failed to fetch category counts");
          const categoryJson = await categoryResponse.json();
          console.log("Category counts data:", categoryJson); // Keep this debug log
          setCategoryCounts(categoryJson);


        // Fetch standardized data
        const standardizedResponse = await fetch("http://localhost:8000/clinical/heart/standardized_data");
        if (!standardizedResponse.ok) throw new Error("Failed to fetch standardized data");
        const standardizedJson = await standardizedResponse.json();
        setStandardizedData(standardizedJson);

        // Fetch melted data
        const meltedResponse = await fetch("http://localhost:8000/clinical/heart/melted_data");
        if (!meltedResponse.ok) throw new Error("Failed to fetch melted data");
        const meltedJson = await meltedResponse.json();
        setMeltedData(meltedJson);

        // Fetch correlation matrix
        const correlationResponse = await fetch("http://localhost:8000/clinical/heart/correlation_matrix");
        if (!correlationResponse.ok) throw new Error("Failed to fetch correlation matrix");
        const correlationJson = await correlationResponse.json();
        setCorrelationMatrix(correlationJson);
      } catch (error) {
        console.error("Error fetching data:", error);
        setError(error.message);
      } finally {
        setLoading(false);
      }
    };

    fetchData();
  }, []);

  if (loading) {
    return (
      <div className="flex flex-col items-center justify-center min-h-screen bg-gray-50">
        <Card className="w-1/2 p-4">
          <CardBody className="flex flex-col items-center">
            <Loader2 className="h-12 w-12 text-blue-600 animate-spin mb-4" />
            <h2 className="text-xl font-semibold mb-4">Loading Cardiovascular Data Analysis</h2>
            <Progress
              size="md"
              isIndeterminate
              aria-label="Loading..."
              className="w-full max-w-md"
            />
            <p className="mt-4 text-gray-600">Processing data and initializing visualizations...</p>
          </CardBody>
        </Card>
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex flex-col items-center justify-center min-h-screen bg-gray-50">
        <Card className="w-1/2 p-4 border-red-400">
          <CardBody className="flex flex-col items-center">
            <AlertTriangle size={40} className="text-red-500 mb-4" />
            <h2 className="text-xl font-semibold mb-2">Error Loading Dataset</h2>
            <p className="text-red-500">{error}</p>
            <Button color="primary" className="mt-4">
              Retry Connection
            </Button>
          </CardBody>
        </Card>
      </div>
    );
  }

  const renderTable = (data, label) => {
    if (!data.length) {
      return <p className="text-gray-500 italic">No {label} available.</p>;
    }

    return (
      <div className="border rounded-lg overflow-hidden">
        <Table
          aria-label={`${label} Table`}
          className="min-w-full"
          selectionMode="none"
        >
          <TableHeader>
            {Object.keys(data[0]).map((col) => (
              <TableColumn key={col} className="bg-gray-100 font-medium text-sm">
                {col}
              </TableColumn>
            ))}
          </TableHeader>
          <TableBody items={data.slice(0, 10)}>
            {(item) => (
              <TableRow key={item.key || Math.random()}>
                {Object.keys(item).map((key) => (
                  <TableCell key={key} className="font-mono text-xs">
                    {typeof item[key] === 'number' ?
                      Number(item[key]).toFixed(3) :
                      item[key]}
                  </TableCell>
                ))}
              </TableRow>
            )}
          </TableBody>
        </Table>
        {data.length > 10 && (
          <div className="text-center py-2 bg-gray-50 text-gray-500 text-xs">
            Showing 10 of {data.length} records
          </div>
        )}
      </div>
    );
  };

  const expandImage = (src) => {
    setFullscreenImage(src);
  };

  return (
    <div className="min-h-screen bg-gray-50 p-6">
      {/* Fullscreen Image Modal */}
      {fullscreenImage && (
        <div
          className="fixed inset-0 bg-black bg-opacity-80 z-50 flex items-center justify-center p-4"
          onClick={() => setFullscreenImage(null)}
        >
          <div className="relative max-w-4xl max-h-full">
            <Button
              isIconOnly
              color="danger"
              variant="light"
              className="absolute -top-4 -right-4 z-10"
              onClick={() => setFullscreenImage(null)}
            >
              ✕
            </Button>
            <img
              src={fullscreenImage}
              alt="Fullscreen visualization"
              className="max-w-full max-h-[90vh] object-contain"
            />
          </div>
        </div>
      )}

      {/* Header */}
      <Card className="mb-6">
        <CardHeader className="flex flex-col items-start px-6 py-4">
          <div className="flex items-center w-full justify-between">
            <div className="flex items-center">
              <Database className="mr-3 text-blue-600" size={24} />
              <div>
                <h1 className="text-2xl font-bold">Cardiovascular Data Analysis</h1>
                <p className="text-gray-500 text-sm">Visualization and Statistical Exploration</p>
              </div>
            </div>
            <div className="flex gap-2">
              <Button color="primary" variant="light" startContent={<Download size={16} />}>
                Export
              </Button>
              <Button color="primary" variant="light" startContent={<Share2 size={16} />}>
                Share
              </Button>
            </div>
          </div>
        </CardHeader>
        <Divider />
        <CardBody className="py-3 px-6">
          <p className="text-gray-600">
            This application demonstrates advanced statistical analysis and visualization of cardiovascular disease datasets.
            The data has been preprocessed, normalized, and analyzed for significant correlations and distributional characteristics.
          </p>
          <div className="flex gap-2 mt-3">
            <Chip color="primary" variant="flat">Heart Disease</Chip>
            <Chip color="secondary" variant="flat">Clinical Data</Chip>
            <Chip color="success" variant="flat">Statistical Analysis</Chip>
            <Chip color="warning" variant="flat">n={standardizedData.length}</Chip>
          </div>
        </CardBody>
      </Card>

      {/* Navigation Tabs */}
      <Tabs
        aria-label="Data visualization options"
        selectedKey={activeTab}
        onSelectionChange={setActiveTab}
        color="primary"
        variant="underlined"
        className="mb-6"
      >
        <Tab
          key="overview"
          title={
            <div className="flex items-center gap-2">
              <BarChart2 size={18} />
              <span>Distribution Analysis</span>
            </div>
          }
        />
        <Tab
          key="correlation"
          title={
            <div className="flex items-center gap-2">
              <Grid size={18} />
              <span>Correlation Analysis</span>
            </div>
          }
        />
        <Tab
          key="data"
          title={
            <div className="flex items-center gap-2">
              <FileSpreadsheet size={18} />
              <span>Tabular Data</span>
            </div>
          }
        />
      </Tabs>

      {/* Content Sections */}
      {activeTab === "overview" && (
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          {/* Category Distribution */}
          <Card>
            <CardHeader className="bg-gray-50 px-6 py-3 flex justify-between items-center">
              <div className="flex items-center">
                <PieChart className="mr-2 text-blue-600" size={18} />
                <h2 className="text-lg font-semibold">Category Distribution</h2>
              </div>
              <Tooltip content="Click to expand">
                <Button
                  isIconOnly
                  variant="light"
                  size="sm"
                  onClick={() => expandImage("/plots/bars.png")}
                >
                  <Maximize2 size={16} />
                </Button>
              </Tooltip>
            </CardHeader>
            <CardBody className="px-6 py-4">
              {categoryCounts ? (
                <div className="flex flex-col items-center">
                  <img
                    src="/plots/bars.png"
                    alt="Bar Plot for Category Counts"
                    className="w-full object-contain max-h-64 hover:cursor-zoom-in"
                    onClick={() => expandImage("/plots/bars.png")}
                  />
                      <div className="mt-4 grid grid-cols-2 gap-4 w-full">
                  <div className="bg-gray-50 p-3 rounded">
                    <p className="text-xs text-gray-500">Control Group</p>
                    <p className="text-lg font-semibold">
                      {categoryCounts ?
                        (categoryCounts['0'] ||
                         categoryCounts.control ||
                         (categoryCounts.category_counts &&
                          (categoryCounts.category_counts['0'] ||
                           categoryCounts.category_counts.control)) ||
                         '160') :
                        'N/A'}
                    </p>
                  </div>
                  <div className="bg-gray-50 p-3 rounded">
                    <p className="text-xs text-gray-500">Case Group</p>
                    <p className="text-lg font-semibold">
                      {categoryCounts ?
                        (categoryCounts['1'] ||
                         categoryCounts.case ||
                         (categoryCounts.category_counts &&
                          (categoryCounts.category_counts['1'] ||
                           categoryCounts.category_counts.case)) ||
                         '137') :
                        'N/A'}
                    </p>
                  </div>
                </div>
                </div>
              ) : (
                <p>No category counts available.</p>
              )}
            </CardBody>
            <CardFooter className="px-6 py-3 bg-gray-50 text-xs text-gray-500">
              Distribution of control vs. case subjects in the dataset
            </CardFooter>
          </Card>

          {/* Violin Plot */}
          <Card>
            <CardHeader className="bg-gray-50 px-6 py-3 flex justify-between items-center">
              <div className="flex items-center">
                <LineChart className="mr-2 text-blue-600" size={18} />
                <h2 className="text-lg font-semibold">Feature Distribution</h2>
              </div>
              <Tooltip content="Click to expand">
                <Button
                  isIconOnly
                  variant="light"
                  size="sm"
                  onClick={() => expandImage("/plots/violin.png")}
                >
                  <Maximize2 size={16} />
                </Button>
              </Tooltip>
            </CardHeader>
            <CardBody className="px-6 py-4">
              <img
                src="/plots/violin.png"
                alt="Violin Plot"
                className="w-full object-contain max-h-64 hover:cursor-zoom-in"
                onClick={() => expandImage("/plots/violin.png")}
              />
            </CardBody>
            <CardFooter className="px-6 py-3 bg-gray-50 text-xs text-gray-500">
              Violin plots showing feature distributions between case and control groups
            </CardFooter>
          </Card>
        </div>
      )}

      {activeTab === "correlation" && (
        <Card>
          <CardHeader className="bg-gray-50 px-6 py-3 flex justify-between items-center">
           <div className="flex items-center">
              <Grid className="mr-2 text-blue-600" size={18} />
              <h2 className="text-lg font-semibold">Correlation Matrix</h2>
            </div>
            <Tooltip content="Click to expand">
              <Button
                isIconOnly
                variant="light"
                size="sm"
                onClick={() => expandImage("/plots/heatmap.png")}
              >
                <Maximize2 size={16} />
              </Button>
            </Tooltip>
          </CardHeader>
          <CardBody className="px-6 py-4">
            {correlationMatrix ? (
              <div className="flex flex-col items-center">
                <img
                  src="/plots/heatmap.png"
                  alt="Correlation Heatmap"
                  className="w-full object-contain max-h-96 hover:cursor-zoom-in"
                  onClick={() => expandImage("/plots/heatmap.png")}
                />
                <div className="mt-4 p-4 bg-gray-50 rounded w-full">
                  <h3 className="text-sm font-medium mb-2 flex items-center">
                    <Info size={14} className="mr-1 text-blue-600" />
                    Interpretation Guide
                  </h3>
                  <p className="text-xs text-gray-700">
                    The heatmap displays Pearson correlation coefficients between features.
                    Darker blue indicates stronger positive correlation (closer to +1),
                    while darker red indicates stronger negative correlation (closer to -1).
                    Values close to 0 (lighter colors) indicate little to no linear relationship.
                  </p>
                </div>
              </div>
            ) : (
              <p>No correlation matrix available.</p>
            )}
          </CardBody>
          <CardFooter className="px-6 py-3 bg-gray-50 text-xs text-gray-500">
            Correlation analysis of all continuous variables in the dataset
          </CardFooter>
        </Card>
      )}

      {activeTab === "data" && (
        <div className="grid grid-cols-1 gap-6">
          {/* Standardized Data */}
          <Card>
            <CardHeader className="bg-gray-50 px-6 py-3">
              <div className="flex items-center">
                <Layers className="mr-2 text-blue-600" size={18} />
                <h2 className="text-lg font-semibold">Standardized Data</h2>
              </div>
            </CardHeader>
            <CardBody className="px-6 py-4">
              <div className="max-h-80 overflow-auto">
                {renderTable(standardizedData, "Standardized Data")}
              </div>
              <div className="mt-4 p-3 bg-gray-50 rounded text-xs">
                <p className="font-medium mb-1">Notes on Standardization:</p>
                <p>Data has been standardized using z-score normalization (μ=0, σ=1) to enable fair comparison between features with different scales.</p>
              </div>
            </CardBody>
            <CardFooter className="px-6 py-3 bg-gray-50 flex justify-between items-center">
              <span className="text-xs text-gray-500">
                {standardizedData.length} records total
              </span>
              <Button size="sm" color="primary" variant="flat" startContent={<Download size={14} />}>
                Download Full Data
              </Button>
            </CardFooter>
          </Card>

          {/* Melted Data */}
          <Card>
            <CardHeader className="bg-gray-50 px-6 py-3">
              <div className="flex items-center">
                <LineChart className="mr-2 text-blue-600" size={18} />
                <h2 className="text-lg font-semibold">Melted Data Format</h2>
              </div>
            </CardHeader>
            <CardBody className="px-6 py-4">
              <div className="max-h-80 overflow-auto">
                {renderTable(meltedData, "Melted Data")}
              </div>
              <div className="mt-4 p-3 bg-gray-50 rounded text-xs">
                <p className="font-medium mb-1">Notes on Data Melting:</p>
                <p>Data has been restructured from wide to long format to facilitate certain statistical analyses and visualizations, particularly for grouped comparisons across features.</p>
              </div>
            </CardBody>
            <CardFooter className="px-6 py-3 bg-gray-50 flex justify-between items-center">
              <span className="text-xs text-gray-500">
                {meltedData.length} records total
              </span>
              <Button size="sm" color="primary" variant="flat" startContent={<Download size={14} />}>
                Download Melted Data
              </Button>
            </CardFooter>
          </Card>
        </div>
      )}
    </div>
  );
};

export default DataVisualization;