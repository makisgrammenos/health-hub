'use client';

import React, { useState } from 'react';
import {
  Card,
  CardHeader,
  CardBody,
  Button,
  Spinner,
  Image,
  // Text,
} from "@nextui-org/react";

const TranscriptomicOverRepresentationAnalysis = () => {
  const [analysisInProgress, setAnalysisInProgress] = useState(false);
  const [analysisOutput, setAnalysisOutput] = useState(null);
  const [executionError, setExecutionError] = useState(null);

  const initiateAnalysisPipeline = async () => {
    setAnalysisInProgress(true);
    setExecutionError(null);
    setAnalysisOutput(null);

    try {
      const apiResponse = await fetch("http://localhost:8000/transcriptomics/ora_analysis", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          "Accept": "application/json"
        }
      });

      if (!apiResponse.ok) {
        const errorPayload = await apiResponse.json();
        throw new Error(errorPayload.detail || "Analysis execution failed: API returned non-200 status");
      }

      const analysisResults = await apiResponse.json();
      setAnalysisOutput(analysisResults);
    } catch (exception) {
      setExecutionError(exception.message);
      console.error("Analysis pipeline error:", exception);
    } finally {
      setAnalysisInProgress(false);
    }
  };

  return (
    <div className="transcriptomics-analysis-container" style={{
      padding: "2rem",
      fontFamily: "Inter, system-ui, sans-serif",
      maxWidth: "1200px",
      margin: "0 auto",
      color: "#1a1a2e"
    }}>
      <header style={{ marginBottom: "2.5rem" }}>
        <h1 style={{
          fontSize: "1.8rem",
          textAlign: "center",
          fontWeight: "600",
          letterSpacing: "-0.025em",
          borderBottom: "1px solid #e2e8f0",
          paddingBottom: "1rem"
        }}>
          Transcriptomic Over-Representation Analysis Pipeline
        </h1>
        <p style={{
          textAlign: "center",
          fontSize: "0.95rem",
          color: "#4a5568",
          maxWidth: "80%",
          margin: "1rem auto 0"
        }}>
          Statistical framework for detecting significantly enriched biological pathways in gene expression data
        </p>
      </header>

      {analysisInProgress ? (
        <div style={{
          textAlign: "center",
          margin: "3rem 0",
          padding: "2rem",
          backgroundColor: "#f7fafc",
          borderRadius: "0.5rem"
        }}>
          <Spinner size="lg" color="primary" />
          <div style={{ marginTop: "1.5rem" }}>
            <p style={{ fontWeight: "500" }}>Processing Transcriptomic Data</p>
            <p style={{
              fontSize: "0.875rem",
              color: "#4a5568",
              marginTop: "0.5rem"
            }}>
              Executing statistical enrichment analysis against reference gene ontology database
            </p>
          </div>
        </div>
      ) : executionError ? (
        <Card style={{
          backgroundColor: "#fff5f5",
          border: "1px solid #feb2b2",
          marginBottom: "2rem"
        }}>
          <CardHeader style={{ backgroundColor: "#fed7d7", paddingTop: "0.75rem", paddingBottom: "0.75rem" }}>
            <h2 style={{ fontSize: "1.1rem", fontWeight: "600", color: "#c53030" }}>
              Analysis Execution Error
            </h2>
          </CardHeader>
          <CardBody>
            <h3 style={{ color: "#c53030" }}>{executionError}</h3>
            <h3 style={{ marginTop: "1rem", fontSize: "0.875rem" }}>
              Check server logs for additional diagnostic information.
            </h3>
          </CardBody>
        </Card>
      ) : analysisOutput ? (
        <div className="analysis-results">
          <Card style={{ marginBottom: "2rem", boxShadow: "0 1px 3px 0 rgba(0, 0, 0, 0.1)" }}>
            <CardHeader style={{
              borderBottom: "1px solid #e2e8f0",
              backgroundColor: "#f7fafc",
              paddingTop: "0.75rem",
              paddingBottom: "0.75rem"
            }}>
              <h2 style={{ fontSize: "1.1rem", fontWeight: "600" }}>Analysis Metadata</h2>
            </CardHeader>
            <CardBody style={{ padding: "1.25rem" }}>
              <div style={{ display: "grid", gridTemplateColumns: "1fr 2fr", gap: "0.5rem" }}>
                <h3 style={{ fontWeight: "500" }}>Status:</h3>
                <h3>{analysisOutput.message}</h3>

                <h3 style={{ fontWeight: "500" }}>Output File:</h3>
                <h3 style={{ fontFamily: "monospace", fontSize: "0.875rem", backgroundColor: "#f1f5f9", padding: "0.25rem 0.5rem", borderRadius: "0.25rem" }}>
                  {analysisOutput.output_file}
                </h3>
              </div>
            </CardBody>
          </Card>

          <div style={{ display: "grid", gridTemplateColumns: "1fr", gap: "2rem" }}>
            <Card style={{ boxShadow: "0 1px 3px 0 rgba(0, 0, 0, 0.1)" }}>
              <CardHeader style={{
                borderBottom: "1px solid #e2e8f0",
                backgroundColor: "#f7fafc",
                paddingTop: "0.75rem",
                paddingBottom: "0.75rem"
              }}>
                <div>
                  <h2 style={{ fontSize: "1.1rem", fontWeight: "600" }}>Leiden Community Detection Clustering</h2>
                  <p style={{ fontSize: "0.8rem", color: "#4a5568", marginTop: "0.25rem" }}>
                    Visualization of transcriptionally similar cell populations
                  </p>
                </div>
              </CardHeader>
              <CardBody style={{ padding: "1rem" }}>
                <div style={{
                  backgroundColor: "#fff",
                  border: "1px solid #e2e8f0",
                  borderRadius: "0.25rem",
                  padding: "1rem"
                }}>
                  <Image
                    src={`data:image/png;base64,${analysisOutput.leiden_plot}`}
                    alt="Leiden algorithm clustering of transcriptionally similar cell populations"
                    objectFit="contain"
                    width="100%"
                  />
                </div>
              </CardBody>
            </Card>

            <Card style={{ boxShadow: "0 1px 3px 0 rgba(0, 0, 0, 0.1)" }}>
              <CardHeader style={{
                borderBottom: "1px solid #e2e8f0",
                backgroundColor: "#f7fafc",
                paddingTop: "0.75rem",
                paddingBottom: "0.75rem"
              }}>
                <div>
                  <h2 style={{ fontSize: "1.1rem", fontWeight: "600" }}>Expression Distribution Analysis</h2>
                  <p style={{ fontSize: "0.8rem", color: "#4a5568", marginTop: "0.25rem" }}>
                    Statistical distribution of gene expression across identified clusters
                  </p>
                </div>
              </CardHeader>
              <CardBody style={{ padding: "1rem" }}>
                <div style={{
                  backgroundColor: "#fff",
                  border: "1px solid #e2e8f0",
                  borderRadius: "0.25rem",
                  padding: "1rem"
                }}>
                  <Image
                    src={`data:image/png;base64,${analysisOutput.violin_plot}`}
                    alt="Violin plot showing distribution of gene expression values across cell clusters"
                    objectFit="contain"
                    width="100%"
                  />
                </div>
              </CardBody>
            </Card>
          </div>
        </div>
      ) : (
        <div style={{
          textAlign: "center",
          margin: "4rem 0",
          padding: "3rem",
          backgroundColor: "#f7fafc",
          borderRadius: "0.5rem",
          border: "1px dashed #cbd5e0"
        }}>
          <h2 style={{
            fontSize: "1.2rem",
            fontWeight: "500",
            marginBottom: "1.5rem"
          }}>
            Ready to Execute Analysis Pipeline
          </h2>
          <p style={{
            fontSize: "0.95rem",
            color: "#4a5568",
            maxWidth: "80%",
            margin: "0 auto 2rem"
          }}>
            The analysis will perform pathway enrichment on the uploaded transcriptomic dataset
          </p>
          <Button
            color="primary"
            size="lg"
            onPress={initiateAnalysisPipeline}
            style={{
              fontWeight: "500",
              padding: "0.5rem 1.5rem"
            }}
          >
            Execute Analysis
          </Button>
        </div>
      )}
    </div>
  );
};

export default TranscriptomicOverRepresentationAnalysis;