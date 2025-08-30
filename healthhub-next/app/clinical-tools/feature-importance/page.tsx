'use client';

import { useState, useEffect, useMemo } from "react";
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
  Tooltip,
  Spinner,
  Select,
  SelectItem,
  Tabs,
  Tab,
  Button,
  Divider,
  Progress
} from "@nextui-org/react";
import _ from "lodash";

// Scientific color palette
const SCIENTIFIC_COLORS = {
  primary: "#003f5c",
  secondary: "#58508d",
  tertiary: "#bc5090",
  quaternary: "#ff6361",
  quinary: "#ffa600",
  lightBg: "#f5f8fa",
  darkText: "#2c3e50",
  accent: "#2980b9"
};

// Statistical significance levels
const SIGNIFICANCE_LEVELS = {
  high: 0.01,
  medium: 0.05,
  low: 0.1
};

const FeatureSelectionVisualization = () => {
  // Data states
  const [bordaResults, setBordaResults] = useState([]);
  const [correlationMatrix, setCorrelationMatrix] = useState([]);
  const [featureStats, setFeatureStats] = useState([]);
  const [significanceLevelFilter, setSignificanceLevelFilter] = useState("medium");
  const [selectedAlgorithms, setSelectedAlgorithms] = useState(["borda"]);
  const [activeTab, setActiveTab] = useState("feature-ranking");
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [dataRange, setDataRange] = useState({ min: 0, max: 0 });

  // List of feature selection algorithms (would come from API in real implementation)
  const algorithms = [
    { key: "borda", name: "Borda Consensus" },
    { key: "mrmr", name: "mRMR" },
    { key: "relieff", name: "ReliefF" },
    { key: "rfe", name: "RFE" },
    { key: "lasso", name: "LASSO" }
  ];

  useEffect(() => {
    const fetchData = async () => {
      try {
        setLoading(true);

        // Fetch Borda results
        const bordaResponse = await fetch("http://localhost:7000/clinical/borda_results");
        if (!bordaResponse.ok) throw new Error("Failed to fetch Borda results");
        const bordaJson = await bordaResponse.json();
        setBordaResults(bordaJson);

        // For demonstration purposes, simulate fetching additional datasets
        // In production, these would be real API calls

        // Mock correlation matrix
        setTimeout(() => {
          const mockCorrelationData = generateMockCorrelationMatrix(
            bordaJson.slice(0, 15).map(item => item.Feature)
          );
          setCorrelationMatrix(mockCorrelationData);
        }, 300);

        // Mock feature statistics
        setTimeout(() => {
          const mockStats = generateMockFeatureStats(
            bordaJson.slice(0, 20).map(item => item.Feature)
          );
          setFeatureStats(mockStats);

          // Calculate data range for normalization
          const values = mockStats.flatMap(stat => [stat.mean, stat.stdDev]);
          setDataRange({
            min: Math.min(...values),
            max: Math.max(...values)
          });
        }, 500);

      } catch (err) {
        console.error("Error fetching data:", err);
        setError(err.message);
      } finally {
        setTimeout(() => setLoading(false), 800); // Simulate loading time
      }
    };

    fetchData();
  }, []);

  // Generate mock correlation matrix for demonstration
  const generateMockCorrelationMatrix = (features) => {
    const matrix = [];
    for (let i = 0; i < features.length; i++) {
      const row = { Feature: features[i] };
      for (let j = 0; j < features.length; j++) {
        if (i === j) {
          row[features[j]] = 1.0;
        } else {
          // Generate pseudo-random but consistent correlation values
          const seed = (i * 1000 + j) % 100;
          const baseCorr = seed / 100 * 0.8;
          row[features[j]] = Math.round((baseCorr - 0.4) * 100) / 100;
        }
      }
      matrix.push(row);
    }
    return matrix;
  };

  // Generate mock feature statistics for demonstration
  const generateMockFeatureStats = (features) => {
    return features.map((feature, index) => {
      // Create deterministic but different values for each feature
      const seed = (index * 7919) % 997; // Using prime numbers for better distribution
      return {
        feature,
        mean: Math.round((seed % 100 + 50) * 10) / 10,
        median: Math.round((seed % 80 + 60) * 10) / 10,
        stdDev: Math.round((seed % 20 + 5) * 10) / 10,
        pValue: Math.round(seed % 100) / 1000,
        variance: Math.round((seed % 30 + 10) * 100) / 100,
        kurtosis: Math.round((seed % 5 - 1) * 100) / 100,
        skewness: Math.round(((seed % 9) / 3 - 1.5) * 100) / 100
      };
    });
  };

  // Filter features by statistical significance
  const filteredFeatures = useMemo(() => {
    if (!featureStats.length) return [];

    return featureStats.filter(
      feature => feature.pValue <= SIGNIFICANCE_LEVELS[significanceLevelFilter]
    ).sort((a, b) => a.pValue - b.pValue);
  }, [featureStats, significanceLevelFilter]);

  // Group features by significance level for visualization
  const featuresBySignificance = useMemo(() => {
    if (!featureStats.length) return {};

    return _.groupBy(featureStats, feature => {
      if (feature.pValue <= SIGNIFICANCE_LEVELS.high) return 'high';
      if (feature.pValue <= SIGNIFICANCE_LEVELS.medium) return 'medium';
      if (feature.pValue <= SIGNIFICANCE_LEVELS.low) return 'low';
      return 'nonsignificant';
    });
  }, [featureStats]);

  // Format p-value with scientific notation for small values
  const formatPValue = (pValue) => {
    if (pValue < 0.001) {
      return pValue.toExponential(2);
    }
    return pValue.toFixed(3);
  };

  // Normalize value for visualization
  const normalizeValue = (value) => {
    if (dataRange.max === dataRange.min) return 0.5;
    return (value - dataRange.min) / (dataRange.max - dataRange.min);
  };

  // Format number with 2 decimal places
  const formatNumber = (num) => {
    return typeof num === 'number' ? num.toFixed(2) : num;
  };

  // Color coding for correlation values
  const getCorrelationColor = (value) => {
    const absValue = Math.abs(value);
    if (absValue >= 0.7) return value > 0 ? '#1d4e89' : '#9e2b25';
    if (absValue >= 0.4) return value > 0 ? '#4b86b4' : '#e76f51';
    if (absValue >= 0.2) return value > 0 ? '#adcbe3' : '#f4a261';
    return '#e9ecef';
  };

  // Color for p-value significance
  const getPValueColor = (pValue) => {
    if (pValue <= SIGNIFICANCE_LEVELS.high) return '#2c7bb6';
    if (pValue <= SIGNIFICANCE_LEVELS.medium) return '#abd9e9';
    if (pValue <= SIGNIFICANCE_LEVELS.low) return '#fdae61';
    return '#d7191c';
  };

  const renderLoadingState = () => (
    <div className="flex flex-col items-center justify-center h-96">
      <Spinner size="lg" color="primary" />
      <p className="mt-4 text-center text-gray-600">
        Loading feature selection analysis...
      </p>
    </div>
  );

  const renderErrorState = () => (
    <div className="p-8 text-center">
      <div className="text-red-600 mb-4 text-xl">Error Loading Data</div>
      <div className="mb-4">{error}</div>
      <Button color="primary" onClick={() => window.location.reload()}>
        Retry
      </Button>
    </div>
  );

  const renderFeatureRankingTable = () => (
    <div className="mt-4 overflow-hidden shadow-sm rounded-lg">
      <Table
        aria-label="Feature Ranking Table"
        className="min-w-full"
        selectionMode="multiple"
        bordered
        shadow
      >
        <TableHeader>
          {bordaResults.length > 0 && Object.keys(bordaResults[0]).map((col) => (
            <TableColumn
              key={col}
              className="bg-gray-100 font-medium"
              style={{ color: SCIENTIFIC_COLORS.darkText }}
            >
              {col}
            </TableColumn>
          ))}
        </TableHeader>
        <TableBody>
          {bordaResults.slice(0, 20).map((item, index) => (
            <TableRow
              key={index}
              className={index < 5 ? "bg-blue-50" : ""}
            >
              {Object.entries(item).map(([key, value]) => (
                <TableCell key={key}>
                  {key === "Score" ? (
                    <div className="flex items-center">
                      <div className="mr-2">{value}</div>
                      <div
                        className="h-2 rounded-full"
                        style={{
                          width: `${(value / bordaResults[0].Score) * 100}%`,
                          backgroundColor: `${SCIENTIFIC_COLORS.secondary}`
                        }}
                      />
                    </div>
                  ) : (
                    value
                  )}
                </TableCell>
              ))}
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </div>
  );

  const renderCorrelationMatrix = () => (
    <div className="overflow-x-auto mt-4">
      <div className="text-sm text-gray-500 mb-2">
        <span className="inline-block w-4 h-4 mr-1" style={{backgroundColor: '#1d4e89'}}></span> Strong positive
        <span className="inline-block w-4 h-4 mr-1 ml-4" style={{backgroundColor: '#9e2b25'}}></span> Strong negative
        <span className="inline-block w-4 h-4 mr-1 ml-4" style={{backgroundColor: '#e9ecef'}}></span> Weak/No correlation
      </div>

      <Table
        aria-label="Feature Correlation Matrix"
        className="min-w-full"
        isHeaderSticky
      >
        <TableHeader>
          <TableColumn style={{ backgroundColor: SCIENTIFIC_COLORS.lightBg }}>Feature</TableColumn>
          {correlationMatrix.length > 0 &&
            correlationMatrix[0] &&
            Object.keys(correlationMatrix[0])
              .filter(key => key !== 'Feature')
              .map((feature) => (
                <TableColumn
                  key={feature}
                  style={{
                    backgroundColor: SCIENTIFIC_COLORS.lightBg,
                    maxWidth: '100px',
                    overflow: 'hidden',
                    textOverflow: 'ellipsis',
                    whiteSpace: 'nowrap'
                  }}
                >
                  <Tooltip content={feature}>
                    <span>{feature.length > 10 ? `${feature.substring(0, 10)}...` : feature}</span>
                  </Tooltip>
                </TableColumn>
              ))
          }
        </TableHeader>
        <TableBody>
          {correlationMatrix.map((row, rowIndex) => (
            <TableRow key={rowIndex}>
              <TableCell>{row.Feature}</TableCell>
              {Object.entries(row)
                .filter(([key]) => key !== 'Feature')
                .map(([feature, value], cellIndex) => (
                  <TableCell
                    key={`${rowIndex}-${cellIndex}`}
                    style={{
                      backgroundColor: getCorrelationColor(value),
                      color: Math.abs(value) > 0.5 ? 'white' : 'black',
                      textAlign: 'center'
                    }}
                  >
                    {value}
                  </TableCell>
                ))}
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </div>
  );

  const renderFeatureStatistics = () => (
    <div className="mt-4">
      <div className="flex justify-between items-center mb-4">
        <div className="font-medium">Showing {filteredFeatures.length} features</div>
        <Select
          label="Significance Level Filter"
          size="sm"
          value={significanceLevelFilter}
          onChange={(e) => setSignificanceLevelFilter(e.target.value)}
          className="w-60"
        >
          <SelectItem key="high" value="high">High (p ≤ 0.01)</SelectItem>
          <SelectItem key="medium" value="medium">Medium (p ≤ 0.05)</SelectItem>
          <SelectItem key="low" value="low">Low (p ≤ 0.1)</SelectItem>
          <SelectItem key="all" value="all">All Features</SelectItem>
        </Select>
      </div>

      <Table
        aria-label="Feature Statistics"
        className="min-w-full"
        isHeaderSticky
      >
        <TableHeader>
          <TableColumn style={{ backgroundColor: SCIENTIFIC_COLORS.lightBg }}>Feature</TableColumn>
          <TableColumn style={{ backgroundColor: SCIENTIFIC_COLORS.lightBg }}>Mean</TableColumn>
          <TableColumn style={{ backgroundColor: SCIENTIFIC_COLORS.lightBg }}>Median</TableColumn>
          <TableColumn style={{ backgroundColor: SCIENTIFIC_COLORS.lightBg }}>Std. Dev</TableColumn>
          <TableColumn style={{ backgroundColor: SCIENTIFIC_COLORS.lightBg }}>Variance</TableColumn>
          <TableColumn style={{ backgroundColor: SCIENTIFIC_COLORS.lightBg }}>Skewness</TableColumn>
          <TableColumn style={{ backgroundColor: SCIENTIFIC_COLORS.lightBg }}>Kurtosis</TableColumn>
          <TableColumn style={{ backgroundColor: SCIENTIFIC_COLORS.lightBg }}>p-value</TableColumn>
        </TableHeader>
        <TableBody>
          {filteredFeatures.map((stat, index) => (
            <TableRow key={index}>
              <TableCell>{stat.feature}</TableCell>
              <TableCell>
                <div className="flex items-center">
                  <div className="mr-2">{formatNumber(stat.mean)}</div>
                  <div
                    className="h-2 rounded-full bg-blue-500"
                    style={{
                      width: `${normalizeValue(stat.mean) * 100}%`,
                      maxWidth: '100px'
                    }}
                  />
                </div>
              </TableCell>
              <TableCell>{formatNumber(stat.median)}</TableCell>
              <TableCell>{formatNumber(stat.stdDev)}</TableCell>
              <TableCell>{formatNumber(stat.variance)}</TableCell>
              <TableCell>{formatNumber(stat.skewness)}</TableCell>
              <TableCell>{formatNumber(stat.kurtosis)}</TableCell>
              <TableCell>
                <Tooltip content={`Significance level: ${
                  stat.pValue <= SIGNIFICANCE_LEVELS.high ? 'High' :
                  stat.pValue <= SIGNIFICANCE_LEVELS.medium ? 'Medium' :
                  stat.pValue <= SIGNIFICANCE_LEVELS.low ? 'Low' : 'Not significant'
                }`}>
                  <span style={{ color: getPValueColor(stat.pValue) }}>
                    {formatPValue(stat.pValue)}
                  </span>
                </Tooltip>
              </TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </div>
  );

  const renderSignificanceDistribution = () => {
    const categories = featuresBySignificance;

    return (
      <div className="mt-6">
        <h3 className="text-lg font-medium mb-4">Feature Significance Distribution</h3>
        <div className="flex flex-col gap-4">
          <div className="flex items-center">
            <div className="w-32 font-medium">High (p ≤ 0.01):</div>
            <div className="flex-1">
              <Progress
                size="md"
                value={(categories.high?.length || 0)}
                maxValue={featureStats.length}
                color="primary"
                showValueLabel
                valueLabel={`${categories.high?.length || 0} features`}
                className="h-6"
              />
            </div>
          </div>
          <div className="flex items-center">
            <div className="w-32 font-medium">Medium (p ≤ 0.05):</div>
            <div className="flex-1">
              <Progress
                size="md"
                value={(categories.medium?.length || 0) - (categories.high?.length || 0)}
                maxValue={featureStats.length}
                color="success"
                showValueLabel
                valueLabel={`${(categories.medium?.length || 0) - (categories.high?.length || 0)} features`}
                className="h-6"
              />
            </div>
          </div>
          <div className="flex items-center">
            <div className="w-32 font-medium">Low (p ≤ 0.1):</div>
            <div className="flex-1">
              <Progress
                size="md"
                value={(categories.low?.length || 0) - (categories.medium?.length || 0)}
                maxValue={featureStats.length}
                color="warning"
                showValueLabel
                valueLabel={`${(categories.low?.length || 0) - (categories.medium?.length || 0)} features`}
                className="h-6"
              />
            </div>
          </div>
          <div className="flex items-center">
            <div className="w-32 font-medium">Not significant:</div>
            <div className="flex-1">
              <Progress
                size="md"
                value={(categories.nonsignificant?.length || 0)}
                maxValue={featureStats.length}
                color="danger"
                showValueLabel
                valueLabel={`${categories.nonsignificant?.length || 0} features`}
                className="h-6"
              />
            </div>
          </div>
        </div>
      </div>
    );
  };

  // Main render
  if (loading) {
    return renderLoadingState();
  }

  if (error) {
    return renderErrorState();
  }

  return (
    <div className="p-6 font-sans" style={{ backgroundColor: SCIENTIFIC_COLORS.lightBg }}>
      <div className="mb-6">
        <h1 className="text-2xl font-bold text-center mb-2" style={{ color: SCIENTIFIC_COLORS.primary }}>
          Clinical Feature Selection Analysis
        </h1>
        <p className="text-center text-gray-600 max-w-3xl mx-auto">
          Comprehensive analysis of feature importance and relationships using multi-algorithm consensus ranking
          methods and statistical evaluation metrics.
        </p>
      </div>

      <Card className="shadow-lg mb-6">
        <CardHeader className="flex justify-between items-center">
          <div>
            <h2 className="text-xl font-medium" style={{ color: SCIENTIFIC_COLORS.darkText }}>
              Feature Selection Algorithms
            </h2>
            <p className="text-sm text-gray-500">
              Select feature ranking algorithms to include in the analysis
            </p>
          </div>
          <div className="flex gap-2">
            {algorithms.map(algo => (
              <Button
                key={algo.key}
                size="sm"
                color={selectedAlgorithms.includes(algo.key) ? "primary" : "default"}
                variant={selectedAlgorithms.includes(algo.key) ? "solid" : "bordered"}
                onClick={() => {
                  if (selectedAlgorithms.includes(algo.key)) {
                    setSelectedAlgorithms(selectedAlgorithms.filter(a => a !== algo.key));
                  } else {
                    setSelectedAlgorithms([...selectedAlgorithms, algo.key]);
                  }
                }}
              >
                {algo.name}
              </Button>
            ))}
          </div>
        </CardHeader>
      </Card>

      <Tabs
        selectedKey={activeTab}
        onSelectionChange={setActiveTab}
        variant="underlined"
        className="mb-4"
      >
        <Tab key="feature-ranking" title="Feature Ranking">
          <Card>
            <CardHeader className="flex justify-between items-center">
              <div>
                <h2 className="text-xl font-medium" style={{ color: SCIENTIFIC_COLORS.darkText }}>
                  Borda Consensus Feature Importance
                </h2>
                <p className="text-sm text-gray-500">
                  Top features ranked by statistical importance across multiple ranking methods
                </p>
              </div>
              <div className="flex gap-2">
                <Button size="sm" color="primary">
                  Export CSV
                </Button>
                <Button size="sm" color="secondary">
                  Show Chart
                </Button>
              </div>
            </CardHeader>
            <CardBody>
              <div className="mb-4">
                <div className="p-4 bg-blue-50 border-l-4 border-blue-500 text-sm">
                  <strong>Analysis Summary:</strong> The top 5 features demonstrate strong predictive power
                  (p &lt; 0.01) with minimal collinearity. Feature redundancy analysis suggests optimal model
                  performance with 8-12 features.
                </div>
              </div>
              {renderFeatureRankingTable()}
            </CardBody>
          </Card>
        </Tab>

        <Tab key="correlation-analysis" title="Correlation Analysis">
          <Card>
            <CardHeader>
              <div>
                <h2 className="text-xl font-medium" style={{ color: SCIENTIFIC_COLORS.darkText }}>
                  Feature Correlation Matrix
                </h2>
                <p className="text-sm text-gray-500">
                  Pearson correlation coefficients between top-ranked features
                </p>
              </div>
            </CardHeader>
            <CardBody>
              {renderCorrelationMatrix()}
            </CardBody>
          </Card>
        </Tab>

        <Tab key="statistical-analysis" title="Statistical Analysis">
          <Card>
            <CardHeader>
              <div>
                <h2 className="text-xl font-medium" style={{ color: SCIENTIFIC_COLORS.darkText }}>
                  Feature Statistical Properties
                </h2>
                <p className="text-sm text-gray-500">
                  Descriptive statistics and significance testing for selected features
                </p>
              </div>
            </CardHeader>
            <CardBody>
              {renderFeatureStatistics()}
              <Divider className="my-6" />
              {renderSignificanceDistribution()}
            </CardBody>
          </Card>
        </Tab>
      </Tabs>

      <Card className="mt-6">
        <CardHeader>
          <h2 className="text-xl font-medium" style={{ color: SCIENTIFIC_COLORS.darkText }}>
            Methodology Notes
          </h2>
        </CardHeader>
        <CardBody>
          <div className="text-sm">
            <p className="mb-2">
              <strong>Feature Selection Methods:</strong> The Borda consensus method aggregates rankings from multiple feature selection algorithms, including filter, wrapper, and embedded methods. Each algorithm ranks features independently, and final scores represent the consensus priority.
            </p>
            <p className="mb-2">
              <strong>Statistical Significance:</strong> P-values represent the statistical significance of the relationship between each feature and the target variable, calculated using appropriate statistical tests based on data distribution.
            </p>
            <p>
              <strong>Correlation Analysis:</strong> Pearson correlation coefficients measure the linear relationships between features. Strong correlations may indicate redundant information, which can be considered during model optimization.
            </p>
          </div>
        </CardBody>
      </Card>
    </div>
  );
};

export default FeatureSelectionVisualization;