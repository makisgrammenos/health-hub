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
} from "@nextui-org/react";

const Home = () => {
  const [rawData, setRawData] = useState([]);
  const [processedData, setProcessedData] = useState([]);
  const [columns, setColumns] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);

  // Fetch dataset from the API
  useEffect(() => {
    const fetchData = async () => {
      try {
        // Fetch raw dataset
        const rawResponse = await fetch("http://127.0.0.1:8000/clinical/heart/dataset/raw");
        if (!rawResponse.ok) {
          throw new Error("Failed to fetch raw dataset");
        }
        const rawJson = await rawResponse.json();
        let fetchedRawData = rawJson.data;

        // Ensure each item has a unique key
        fetchedRawData = fetchedRawData.map((item, index) => ({
          ...item,
          key: item.id || index, // Replace 'id' with your unique identifier if available
        }));
        setRawData(fetchedRawData);

        // Transform columns to match NextUI's expected structure
        const transformedColumns = rawJson.columns.map((col) => ({
          key: col,
          label: col.toUpperCase(), // Customize label as needed
        }));
        setColumns(transformedColumns);

        // Fetch processed dataset
        const processedResponse = await fetch("http://127.0.0.1:8000/clinical/heart/dataset/processed");
        if (!processedResponse.ok) {
          throw new Error("Failed to fetch processed dataset");
        }
        const processedJson = await processedResponse.json();
        let fetchedProcessedData = processedJson.data;

        // Ensure each item has a unique key
        fetchedProcessedData = fetchedProcessedData.map((item, index) => ({
          ...item,
          key: item.id || index, // Replace 'id' with your unique identifier if available
        }));
        setProcessedData(fetchedProcessedData);
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
    return <div>Loading datasets...</div>;
  }

  if (error) {
    return <div>Error: {error}</div>;
  }

  // Define styles for the scrollable container
  const scrollableContainerStyle = {
    maxHeight: "400px", // Adjust the height as needed
    overflowY: "auto",
    overflowX: "auto",
    // Optional: Add border and padding for better aesthetics
    border: "1px solid #eaeaea",
    borderRadius: "8px",
    padding: "1rem",
  };

  // Create a map for raw data for quick lookup
  const rawDataMap = new Map();
  rawData.forEach((item) => {
    rawDataMap.set(item.key, item);
  });

  // Legend for highlights
  const legendStyle = {
    display: "flex",
    alignItems: "center",
    marginBottom: "1rem",
  };

  const legendItemStyle = {
    display: "flex",
    alignItems: "center",
    marginRight: "1rem",
  };

  const colorBoxStyle = (color) => ({
    width: "16px",
    height: "16px",
    backgroundColor: color,
    marginRight: "0.5rem",
    borderRadius: "4px",
  });

  return (
    <div style={{ padding: "2rem" }}>
      <h1>Dataset Demonstration</h1>
      <p style={{ fontSize: "18px", marginBottom: "1rem" }}>
        This application demonstrates the preprocessing done on the Heart Disease
        dataset. Below, you can see the raw dataset and its processed version.
      </p>

      {/* Legend for Highlights */}
      <div style={legendStyle}>
        <div style={legendItemStyle}>
          <div style={colorBoxStyle("rgba(255, 235, 59, 0.4)")}></div>
          <span>Changed Value</span>
        </div>
        <div style={legendItemStyle}>
          <div style={colorBoxStyle("rgba(76, 175, 80, 0.4)")}></div>
          <span>New Value</span>
        </div>
      </div>

      {/* Raw Dataset Table */}
      <h2>Raw Dataset</h2>
      <div style={scrollableContainerStyle}>
        <Table
          aria-label="Raw Dataset Table"
          css={{
            minWidth: "100%", // Ensures table takes full width of container
          }}
          selectionMode="none" // Optional: Disable row selection if not needed
        >
          <TableHeader columns={columns}>
            {(column) => (
              <TableColumn
                key={column.key}
                css={{
                  position: "sticky",
                  top: 0,
                  backgroundColor: "$accents0", // Adjust based on your theme
                  zIndex: 1,
                }}
              >
                {column.label}
              </TableColumn>
            )}
          </TableHeader>
          <TableBody items={rawData}>
            {(item) => (
              <TableRow key={item.key}>
                {(columnKey) => (
                  <TableCell>
                    {getKeyValue(item, columnKey)}
                  </TableCell>
                )}
              </TableRow>
            )}
          </TableBody>
        </Table>
      </div>

      {/* Processed Dataset Table */}
      <h2 style={{ marginTop: "2rem" }}>Processed Dataset</h2>
      <div style={scrollableContainerStyle}>
        <Table
          aria-label="Processed Dataset Table"
          css={{
            minWidth: "100%", // Ensures table takes full width of container
          }}
          selectionMode="none" // Optional: Disable row selection if not needed
        >
          <TableHeader columns={columns}>
            {(column) => (
              <TableColumn
                key={column.key}
                css={{
                  position: "sticky",
                  top: 0,
                  backgroundColor: "$accents0", // Adjust based on your theme
                  zIndex: 1,
                }}
              >
                {column.label}
              </TableColumn>
            )}
          </TableHeader>
          <TableBody items={processedData}>
            {(item) => (
              <TableRow key={item.key}>
                {(columnKey) => {
                  const rawItem = rawDataMap.get(item.key);
                  const rawValue = rawItem ? rawItem[columnKey] : null;
                  const processedValue = getKeyValue(item, columnKey);
                  const isNew = !rawItem; // If there's no corresponding raw item
                  const isDifferent = rawValue !== processedValue;
                  let backgroundColor = "inherit";
                  let tooltipContent = "";

                  if (isNew) {
                    backgroundColor = "rgba(76, 175, 80, 0.4)"; // Light green for new values
                    tooltipContent = "New Value";
                  } else if (isDifferent) {
                    backgroundColor = "rgba(255, 235, 59, 0.4)"; // Light yellow for changed values
                    tooltipContent = "Changed Value";
                  }

                  return (
                    <TableCell
                      css={{
                        backgroundColor,
                        transition: "background-color 0.3s ease",
                      }}
                    >
                      {isDifferent || isNew ? (
                        <Tooltip content={tooltipContent}>
                          <span>{processedValue}</span>
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
    </div>
  );
};

export default Home;
