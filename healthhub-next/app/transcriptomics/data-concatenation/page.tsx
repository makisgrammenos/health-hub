'use client';

import { useState, useEffect } from "react";
import {
  Card,
  CardHeader,
  CardBody,
  Button,
  Spinner,
  Table,
  TableHeader,
  TableBody,
  TableRow,
  TableCell,
  TableColumn,
} from "@nextui-org/react";

const PreprocessingDemo = () => {
  const [loading, setLoading] = useState(false);
  const [summaries, setSummaries] = useState([]);
  const [samples, setSamples] = useState([]);
  const [error, setError] = useState(null);

  const fetchSummaries = async () => {
    setLoading(true);
    try {
      const response = await fetch("http://localhost:8000/transcriptomics/preprocessing/summary");
      if (!response.ok) throw new Error("Failed to fetch summary");
      const data = await response.json();
      setSummaries(data.summaries);
    } catch (err) {
      console.error("Error fetching summaries:", err);
      setError(err.message);
    } finally {
      setLoading(false);
    }
  };

  const fetchSamples = async () => {
    setLoading(true);
    try {
      const response = await fetch("http://localhost:8000/transcriptomics/preprocessing/sample?rows=5&cols=5");
      if (!response.ok) throw new Error("Failed to fetch samples");
      const data = await response.json();
      setSamples(data.samples);
    } catch (err) {
      console.error("Error fetching samples:", err);
      setError(err.message);
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    fetchSummaries();
  }, []);

  const renderTable = (data, label) => (
    <Table>
      <TableHeader>
        {data.columns.map((col) => (
          <TableColumn key={col}>{col}</TableColumn>
        ))}
      </TableHeader>
      <TableBody>
        {data.data.map((row, rowIndex) => (
          <TableRow key={rowIndex}>
            {row.map((cell, cellIndex) => (
              <TableCell key={cellIndex}>{cell}</TableCell>
            ))}
          </TableRow>
        ))}
      </TableBody>
    </Table>
  );

  return (
    <div style={{ padding: "2rem", fontFamily: "Arial, sans-serif", lineHeight: "1.5" }}>
      <h1 style={{ textAlign: "center", fontSize: "2rem", marginBottom: "2rem" }}>
        Data Preprocessing Demonstration
      </h1>
      <p style={{ textAlign: "center", fontSize: "1rem", marginBottom: "2rem" }}>
        This demonstration highlights preprocessing of raw datasets into a filtered and optimized format.
      </p>

      {loading ? (
        <div style={{ textAlign: "center", margin: "2rem 0" }}>
          <Spinner size="lg" />
          <p style={{ marginTop: "1rem", fontSize: "1.2rem" }}>Loading, please wait...</p>
        </div>
      ) : error ? (
        <div style={{ textAlign: "center", color: "red" }}>
          <p>Error: {error}</p>
        </div>
      ) : (
        <>
          <Card>
            <CardHeader>
              <h2>Summaries</h2>
            </CardHeader>
            <CardBody>
              {summaries.map((summary, index) => (
                <div key={index} style={{ marginBottom: "2rem" }}>
                  <h3>File Pair: {summary.file_pair.raw} → {summary.file_pair.processed}</h3>
                  <p><strong>Raw Data:</strong></p>
                  <p>Cells: {summary.raw_summary.n_cells}, Genes: {summary.raw_summary.n_genes}, Sparsity: {summary.raw_summary.sparsity?.toFixed(2)}</p>
                  <p><strong>Processed Data:</strong></p>
                  <p>Cells: {summary.processed_summary.n_cells}, Genes: {summary.processed_summary.n_genes}, Sparsity: {summary.processed_summary.sparsity?.toFixed(2)}</p>
                </div>
              ))}
            </CardBody>
          </Card>

          <Button style={{ margin: "2rem auto", display: "block" }} onPress={fetchSamples}>
            Fetch Sample Data
          </Button>

          {samples.length > 0 && (
            <Card>
              <CardHeader>
                <h2>Sample Data</h2>
              </CardHeader>
              <CardBody>
                {samples.map((sample, index) => (
                  <div key={index} style={{ marginBottom: "2rem" }}>
                    <h3>File Pair: {sample.file_pair.raw} → {sample.file_pair.processed}</h3>
                    <p><strong>Raw Sample:</strong></p>
                    {renderTable(sample.raw_sample, "Raw Data Sample")}
                    <p><strong>Processed Sample:</strong></p>
                    {renderTable(sample.processed_sample, "Processed Data Sample")}
                  </div>
                ))}
              </CardBody>
            </Card>
          )}
        </>
      )}
    </div>
  );
};

export default PreprocessingDemo;
