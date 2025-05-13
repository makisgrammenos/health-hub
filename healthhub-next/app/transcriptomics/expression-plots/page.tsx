"use client";

import { useState } from "react";
import {
  Button,
  Card,
  CardHeader,
  CardBody,
  Spinner,
  Image,
  Divider,
} from "@nextui-org/react";

export default function UMAPVisualization() {
  const [loading, setLoading] = useState(false);
  const [imageBase64, setImageBase64] = useState("");

  const fetchUMAPImage = async () => {
    setLoading(true);
    try {
      const response = await fetch("http://localhost:8000/transcriptomics/umap-expression");
      if (!response.ok) {
        throw new Error("Failed to fetch the UMAP image.");
      }
      const data = await response.json();
      setImageBase64(data.image_base64);
    } catch (error) {
      console.error("Error fetching UMAP image:", error);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div style={{ padding: "2rem", fontFamily: "Arial, sans-serif" }}>
      <h1 style={{ fontSize: "2.5rem", marginBottom: "1.5rem", textAlign: "center" }}>
        UMAP Gene Expression Visualization
      </h1>
      <p style={{ fontSize: "1.2rem", textAlign: "center", marginBottom: "2rem", color: "#555" }}>
        This demonstration showcases a UMAP (Uniform Manifold Approximation and Projection) visualization 
        of selected gene expressions across a dataset. UMAP is a dimensionality reduction technique widely 
        used in biological and genomic studies to visualize patterns in high-dimensional data.
      </p>
      <Divider style={{ margin: "1rem 0" }} />

      {loading ? (
        <div style={{ textAlign: "center", margin: "2rem 0" }}>
          <Spinner size="lg" />
          <p style={{ fontSize: "1.1rem", color: "#555" }}>
            Generating the visualization. This may take a few moments...
          </p>
        </div>
      ) : imageBase64 ? (
        <Card style={{ maxWidth: "800px", margin: "2rem auto" }}>
          <CardHeader>
            <h2 style={{ fontSize: "1.8rem", textAlign: "center", margin: "1rem 0" }}>
              UMAP Visualization for Selected Genes
            </h2>
          </CardHeader>
          <CardBody>
            <p style={{ textAlign: "center", marginBottom: "1rem", fontSize: "1rem", color: "#555" }}>
              Below is the UMAP visualization plot highlighting the expression of selected genes. 
              Each point represents a cell, and the colors denote the expression levels of the genes 
              across the dataset.
            </p>
            <Image
              src={`data:image/png;base64,${imageBase64}`}
              alt="UMAP Expression Plot"
              style={{
                width: "100%",
                borderRadius: "8px",
                boxShadow: "0 4px 12px rgba(0, 0, 0, 0.1)",
              }}
            />
          </CardBody>
        </Card>
      ) : (
        <div style={{ textAlign: "center" }}>
          <Card style={{ maxWidth: "600px", margin: "0 auto", padding: "1rem" }}>
            <CardHeader>
              <h2 style={{ fontSize: "1.5rem", marginBottom: "0.5rem" }}>What You Will See</h2>
            </CardHeader>
            <CardBody>
              <p style={{ fontSize: "1rem", color: "#555" }}>
                Click the button below to generate the UMAP plot for a set of selected genes. This 
                will demonstrate how gene expression patterns vary across the dataset using the 
                UMAP dimensionality reduction method.
              </p>
              <p style={{ fontSize: "1rem", color: "#555" }}>
                The plot will be dynamically generated and displayed below for review and exploration.
              </p>
            </CardBody>
          </Card>
          <Button
            onClick={fetchUMAPImage}
            style={{
              marginTop: "2rem",
              padding: "0.75rem 2rem",
              backgroundColor: "#0070f3",
              color: "white",
              borderRadius: "5px",
              fontSize: "1rem",
              cursor: "pointer",
            }}
          >
            Generate UMAP Visualization
          </Button>
        </div>
      )}
    </div>
  );
}
