'use client';

import { useState } from "react";
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

const IntegrationDemo = () => {
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState(null);

  const handleIntegration = async () => {
    setLoading(true);
    try {
      const response = await fetch("http://localhost:8000/transcriptomics/integrate", {
        method: "POST",
      });
      if (!response.ok) throw new Error("Integration failed.");
      const data = await response.json();
      setResult(data);
    } catch (error) {
      console.error("Error:", error);
      alert("Integration failed. Please try again.");
    } finally {
      setLoading(false);
    }
  };

  return (
    <div style={{ padding: "2rem", fontFamily: "Arial, sans-serif", lineHeight: "1.5" }}>
      <h1 style={{ textAlign: "center", fontSize: "2rem", marginBottom: "2rem" }}>
        Data Integration Demo
      </h1>

      {loading ? (
        <div style={{ textAlign: "center", margin: "2rem 0" }}>
          <Spinner size="lg" />
          <p>Processing data, please wait...</p>
        </div>
      ) : result ? (
        <Card>
          <CardHeader>
            <h2>Integration Results</h2>
          </CardHeader>
          <CardBody>
            <p>
              <strong>Output File:</strong> {result.output_file}
            </p>
            <p>
              <strong>Batches:</strong>
            </p>
            <Table>
              <TableHeader>
                <TableColumn>Batch</TableColumn>
                <TableColumn>Cell Count</TableColumn>
              </TableHeader>
              <TableBody>
                {Object.entries(result.batches).map(([batch, count], index) => (
                  <TableRow key={index}>
                    <TableCell>{batch}</TableCell>
                    <TableCell>{count}</TableCell>
                  </TableRow>
                ))}
              </TableBody>
            </Table>
            <p style={{ marginTop: "1rem" }}>
              Integration completed successfully. You can download the integrated file from the
              provided path.
            </p>
          </CardBody>
        </Card>
      ) : (
        <Button
          style={{ display: "block", margin: "2rem auto" }}
          onPress={handleIntegration}
        >
          Start Integration
        </Button>
      )}
    </div>
  );
};

export default IntegrationDemo;
