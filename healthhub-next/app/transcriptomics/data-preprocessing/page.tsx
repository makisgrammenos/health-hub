'use client';

import { useState, useEffect } from "react";
import {
  Card,
  CardHeader,
  CardBody,
  Table,
  TableHeader,
  TableColumn,
  TableBody,
  TableRow,
  TableCell,
  Image,
} from "@nextui-org/react";

const DataPreprocessingDemo = () => {
  const [summaries, setSummaries] = useState([]);
  const [samples, setSamples] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);

  useEffect(() => {
    const fetchData = async () => {
      try {
        // Fetch summaries
        const summaryResponse = await fetch("http://localhost:8000/transcriptomics/preprocessing/summary");
        if (!summaryResponse.ok) throw new Error("Failed to fetch summaries");
        const summaryJson = await summaryResponse.json();
        setSummaries(summaryJson.summaries);

        // Fetch samples
        const sampleResponse = await fetch("http://localhost:8000/transcriptomics/preprocessing/sample?rows=10&cols=10");
        if (!sampleResponse.ok) throw new Error("Failed to fetch samples");
        const sampleJson = await sampleResponse.json();
        setSamples(sampleJson.samples);
      } catch (err) {
        console.error("Error fetching data:", err);
        setError(err.message);
      } finally {
        setLoading(false);
      }
    };

    fetchData();
  }, []);

  if (loading) {
    return <div>Loading data...</div>;
  }

  if (error) {
    return <div>Error: {error}</div>;
  }

  return (
    <div style={{ padding: "2rem", fontFamily: "Arial, sans-serif", lineHeight: "1.5" }}>
      <h1 style={{ textAlign: "center", fontSize: "2rem", marginBottom: "2rem" }}>
        Data Preprocessing Demonstration
      </h1>

      {summaries.map((summary, index) => (
        <Card key={index} style={{ marginBottom: "2rem" }}>
          <CardHeader>
            <h2 style={{ margin: 0 }}>File Pair: {summary.file_pair.raw} &rarr; {summary.file_pair.processed}</h2>
          </CardHeader>
          <CardBody>
            <div style={{ display: "flex", flexWrap: "wrap", gap: "1rem" }}>
              <div style={{ flex: 1 }}>
                <h3>Original Data</h3>
                <p>Number of Cells: {summary.raw_summary.n_cells}</p>
                <p>Number of Genes: {summary.raw_summary.n_genes}</p>
                <p>Sparsity: {summary.raw_summary.sparsity ? summary.raw_summary.sparsity.toFixed(2) : "N/A"}</p>
              </div>
              <div style={{ flex: 1 }}>
                <h3>Processed Data</h3>
                <p>Number of Cells: {summary.processed_summary.n_cells}</p>
                <p>Number of Genes: {summary.processed_summary.n_genes}</p>
                <p>Sparsity: {summary.processed_summary.sparsity ? summary.processed_summary.sparsity.toFixed(2) : "N/A"}</p>
              </div>
            </div>
          </CardBody>
        </Card>
      ))}

      {samples.map((sample, index) => (
        <Card key={index} style={{ marginBottom: "2rem" }}>
          <CardHeader>
            <h2 style={{ margin: 0 }}>Sample View: {sample.file_pair.raw} &rarr; {sample.file_pair.processed}</h2>
          </CardHeader>
          <CardBody>
            <h3>Original Data Sample</h3>
            <Table
              aria-label="Original Data Sample"
              style={{ marginBottom: "1rem" }}
            >
              <TableHeader>
                {sample.raw_sample.columns.map((col) => (
                  <TableColumn key={col}>{col}</TableColumn>
                ))}
              </TableHeader>
              <TableBody>
                {sample.raw_sample.data.map((row, rowIndex) => (
                  <TableRow key={rowIndex}>
                    {row.map((cell, cellIndex) => (
                      <TableCell key={cellIndex}>{cell}</TableCell>
                    ))}
                  </TableRow>
                ))}
              </TableBody>
            </Table>

            <h3>Processed Data Sample</h3>
            <Table aria-label="Processed Data Sample">
              <TableHeader>
                {sample.processed_sample.columns.map((col) => (
                  <TableColumn key={col}>{col}</TableColumn>
                ))}
              </TableHeader>
              <TableBody>
                {sample.processed_sample.data.map((row, rowIndex) => (
                  <TableRow key={rowIndex}>
                    {row.map((cell, cellIndex) => (
                      <TableCell key={cellIndex}>{cell}</TableCell>
                    ))}
                  </TableRow>
                ))}
              </TableBody>
            </Table>
          </CardBody>
        </Card>
      ))}
    </div>
  );
};

export default DataPreprocessingDemo;
