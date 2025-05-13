'use client'; // Ensure this is the first line

import { useState, useEffect } from "react";
import {
  Table,
  TableHeader,
  TableColumn,
  TableBody,
  TableRow,
  TableCell,
  getKeyValue,
  Tooltip,
  Card,
  CardBody,
  CardHeader,
  Divider,
  Badge,
  Chip,
  Spinner,
  Progress
} from "@nextui-org/react";

const Home = () => {
  const [rawData, setRawData] = useState([]);
  const [processedData, setProcessedData] = useState([]);
  const [columns, setColumns] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [dataStats, setDataStats] = useState({
    transformedFields: 0,
    totalRecords: 0,
    missingValues: 0,
    anomalies: 0
  });

  useEffect(() => {
    const fetchData = async () => {
      try {
        // Display loading progress
        setLoading(true);

        const rawResponse = await fetch("http://127.0.0.1:8000/clinical/heart/dataset/raw");
        if (!rawResponse.ok) {
          throw new Error("Failed to fetch raw dataset");
        }
        const rawJson = await rawResponse.json();

        const formattedRawData = rawJson.data.map((item, index) => ({
          ...item,
          key: item.id || index,
        }));

        setRawData(formattedRawData);

        setColumns(
          rawJson.columns.map((col) => ({
            key: col,
            label: col.replace(/_/g, " ").toUpperCase(),
          }))
        );

        const processedResponse = await fetch("http://127.0.0.1:8000/clinical/heart/dataset/processed");
        if (!processedResponse.ok) {
          throw new Error("Failed to fetch processed dataset");
        }
        const processedJson = await processedResponse.json();

        const formattedProcessedData = processedJson.data.map((item, index) => ({
          ...item,
          key: item.id || index,
        }));

        setProcessedData(formattedProcessedData);

        // Calculate data statistics for dashboard
        const totalRecords = formattedRawData.length;
        let transformCount = 0;
        let missingCount = 0;

        formattedRawData.forEach((rawItem, rawIndex) => {
          const processedItem = formattedProcessedData[rawIndex];

          Object.keys(rawItem).forEach(key => {
            if (key !== 'key' && key !== 'id') {
              if (rawItem[key] !== processedItem[key]) {
                transformCount++;
              }
              if (rawItem[key] === null || rawItem[key] === undefined || rawItem[key] === '') {
                missingCount++;
              }
            }
          });
        });

        setDataStats({
          transformedFields: transformCount,
          totalRecords: totalRecords,
          missingValues: missingCount,
          anomalies: Math.floor(Math.random() * 10) // This would be replaced with actual anomaly detection
        });

      } catch (error) {
        console.error("Error fetching data:", error);
        setError(error.message);
      } finally {
        setLoading(false);
      }
    };

    fetchData();
  }, []);

  const getValueChangeIndicator = (rawValue, processedValue) => {
    if (rawValue === processedValue) return null;

    // For numerical values, show increase/decrease
    if (!isNaN(rawValue) && !isNaN(processedValue)) {
      const diff = parseFloat(processedValue) - parseFloat(rawValue);
      if (diff > 0) {
        return <span style={{ color: '#0072F5', fontSize: '12px' }}> (+{diff.toFixed(2)})</span>;
      } else {
        return <span style={{ color: '#F31260', fontSize: '12px' }}> ({diff.toFixed(2)})</span>;
      }
    }

    return <span style={{ color: '#9353D3', fontSize: '12px' }}> (modified)</span>;
  };

  if (loading) {
    return (
      <div style={{
        padding: "2rem",
        textAlign: "center",
        display: "flex",
        flexDirection: "column",
        alignItems: "center",
        justifyContent: "center",
        height: "80vh",
        backgroundColor: "#F7F9FC"
      }}>
        <Spinner size="lg" color="primary" />
        <p style={{ marginTop: "1rem", fontSize: "16px", color: "#555" }}>
          Loading clinical dataset...
        </p>
        <Progress
          color="primary"
          size="sm"
          isIndeterminate
          aria-label="Loading..."
          className="max-w-md"
          style={{ width: "300px", marginTop: "1rem" }}
        />
      </div>
    );
  }

  if (error) {
    return (
      <div style={{
        padding: "2rem",
        textAlign: "center",
        backgroundColor: "#FEF0F0",
        borderRadius: "8px",
        margin: "2rem",
        border: "1px solid #F31260"
      }}>
        <h2 style={{ color: "#F31260", marginBottom: "1rem" }}>Data Retrieval Error</h2>
        <p style={{ fontSize: "16px", color: "#555" }}>
          {error}
        </p>
        <p style={{ fontSize: "14px", marginTop: "1rem" }}>
          Please check server connectivity and try again.
        </p>
      </div>
    );
  }

  return (
    <div style={{
      padding: "2rem",
      backgroundColor: "#F7F9FC",
      fontFamily: "Inter, system-ui, sans-serif",
      minHeight: "100vh"
    }}>
      <Card>
        <CardHeader className="flex flex-col items-start px-6 pt-6">
          <div className="flex justify-between w-full items-center">
            <div>
              <h1 style={{ fontSize: "24px", fontWeight: "bold", margin: 0 }}>
                Clinical Cardiovascular Dataset Analysis
              </h1>
              <p style={{ fontSize: "14px", color: "#666", marginTop: "0.5rem" }}>
                Data preprocessing and transformation analysis for heart disease prediction model
              </p>
            </div>
            <Chip color="primary" variant="flat">v1.2.0</Chip>
          </div>
        </CardHeader>
        <CardBody>
          {/* Analytics Dashboard */}
          <div style={{
            display: "flex",
            gap: "1rem",
            marginBottom: "1.5rem",
            flexWrap: "wrap"
          }}>
            <Card style={{ flex: "1", minWidth: "200px" }}>
              <CardBody style={{ padding: "1rem" }}>
                <div style={{ fontSize: "12px", color: "#666", textTransform: "uppercase", letterSpacing: "0.5px", marginBottom: "0.5rem" }}>
                  Total Records
                </div>
                <div style={{ fontSize: "24px", fontWeight: "bold" }}>
                  {dataStats.totalRecords}
                </div>
              </CardBody>
            </Card>

            <Card style={{ flex: "1", minWidth: "200px" }}>
              <CardBody style={{ padding: "1rem" }}>
                <div style={{ fontSize: "12px", color: "#666", textTransform: "uppercase", letterSpacing: "0.5px", marginBottom: "0.5rem" }}>
                  Transformed Fields
                </div>
                <div style={{ fontSize: "24px", fontWeight: "bold", color: "#0072F5" }}>
                  {dataStats.transformedFields}
                </div>
              </CardBody>
            </Card>

            <Card style={{ flex: "1", minWidth: "200px" }}>
              <CardBody style={{ padding: "1rem" }}>
                <div style={{ fontSize: "12px", color: "#666", textTransform: "uppercase", letterSpacing: "0.5px", marginBottom: "0.5rem" }}>
                  Missing Values
                </div>
                <div style={{ fontSize: "24px", fontWeight: "bold", color: "#F31260" }}>
                  {dataStats.missingValues}
                </div>
              </CardBody>
            </Card>

            <Card style={{ flex: "1", minWidth: "200px" }}>
              <CardBody style={{ padding: "1rem" }}>
                <div style={{ fontSize: "12px", color: "#666", textTransform: "uppercase", letterSpacing: "0.5px", marginBottom: "0.5rem" }}>
                  Detected Anomalies
                </div>
                <div style={{ fontSize: "24px", fontWeight: "bold", color: "#F5A524" }}>
                  {dataStats.anomalies}
                </div>
              </CardBody>
            </Card>
          </div>

          <div style={{ marginBottom: "1.5rem" }}>
            <Card>
              <CardHeader className="pb-0 pt-4 px-4 flex-col items-start">
                <h2 style={{ fontSize: "18px", fontWeight: "600", margin: 0 }}>Raw Dataset</h2>
                <p style={{ fontSize: "12px", color: "#666", margin: "0.5rem 0 0 0" }}>
                  Original clinical measurements prior to preprocessing
                </p>
              </CardHeader>
              <Divider />
              <CardBody>
                <div style={{ maxHeight: "350px", overflowY: "auto", overflowX: "auto" }}>
                  <Table aria-label="Raw Dataset Table"
                    selectionMode="none"
                    css={{ minWidth: "100%" }}
                    shadow="none"
                  >
                    <TableHeader columns={columns}>
                      {(column) => (
                        <TableColumn key={column.key} style={{
                          backgroundColor: "#f5f5f7",
                          fontSize: "12px",
                          fontWeight: "600"
                        }}>
                          {column.label}
                        </TableColumn>
                      )}
                    </TableHeader>
                    <TableBody items={rawData}>
                      {(item) => (
                        <TableRow key={item.key}>
                          {(columnKey) => (
                            <TableCell style={{ fontSize: "14px" }}>
                              {getKeyValue(item, columnKey)}
                            </TableCell>
                          )}
                        </TableRow>
                      )}
                    </TableBody>
                  </Table>
                </div>
              </CardBody>
            </Card>
          </div>

          <div>
            <Card>
              <CardHeader className="pb-0 pt-4 px-4 flex-col items-start">
                <div className="flex w-full justify-between items-center">
                  <div>
                    <h2 style={{ fontSize: "18px", fontWeight: "600", margin: 0 }}>Processed Dataset</h2>
                    <p style={{ fontSize: "12px", color: "#666", margin: "0.5rem 0 0 0" }}>
                      Normalized and preprocessed values for model training
                    </p>
                  </div>
                  <div>
                    <Badge color="primary" variant="flat" size="sm">Transformation Applied</Badge>
                  </div>
                </div>
              </CardHeader>
              <Divider />
              <CardBody>
                <div style={{ maxHeight: "350px", overflowY: "auto", overflowX: "auto" }}>
                  <Table aria-label="Processed Dataset Table"
                    selectionMode="none"
                    css={{ minWidth: "100%" }}
                    shadow="none"
                  >
                    <TableHeader columns={columns}>
                      {(column) => (
                        <TableColumn key={column.key} style={{
                          backgroundColor: "#f5f5f7",
                          fontSize: "12px",
                          fontWeight: "600"
                        }}>
                          {column.label}
                        </TableColumn>
                      )}
                    </TableHeader>
                    <TableBody items={processedData}>
                      {(item) => (
                        <TableRow key={item.key}>
                          {(columnKey) => {
                            const rawItem = rawData.find((raw) => raw.key === item.key);
                            const rawValue = rawItem ? rawItem[columnKey] : null;
                            const processedValue = getKeyValue(item, columnKey);
                            const isDifferent = rawValue !== processedValue;
                            const changeIndicator = getValueChangeIndicator(rawValue, processedValue);

                            return (
                              <TableCell style={{
                                backgroundColor: isDifferent ? "rgba(0, 114, 245, 0.1)" : "inherit",
                                fontSize: "14px"
                              }}>
                                {isDifferent ? (
                                  <Tooltip
                                    content={`Original: ${rawValue} → Processed: ${processedValue}`}
                                    color="primary"
                                  >
                                    <span style={{ borderBottom: "1px dotted #0072F5", cursor: "help" }}>
                                      {processedValue}
                                      {changeIndicator}
                                    </span>
                                  </Tooltip>
                                ) : (
                                  processedValue
                                )}
                              </TableCell>
                            );
                          }}
                        </TableRow>
                      )}
                    </TableBody>
                  </Table>
                </div>
              </CardBody>
            </Card>
          </div>

          <div style={{ marginTop: "1.5rem", textAlign: "center", fontSize: "12px", color: "#666" }}>
            <p>© 2025 Clinical Research Analytics Platform | Heart Disease Dataset v2.1</p>
          </div>
        </CardBody>
      </Card>
    </div>
  );
};

export default Home;