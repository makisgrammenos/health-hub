'use client';
import { useState } from "react";
import axios from "axios";
import { Button, Input, Textarea, Card, Loading, Checkbox, Dropdown } from "@nextui-org/react";

export default function DatasetProcessor() {
  const [file, setFile] = useState(null);
  const [columns, setColumns] = useState([]);
  const [isFirstLineHeader, setIsFirstLineHeader] = useState(true);
  const [dropMissing, setDropMissing] = useState(false);
  const [fillColumn, setFillColumn] = useState(null);
  const [fillValue, setFillValue] = useState("");
  const [categoricalColumn, setCategoricalColumn] = useState(null);
  const [threshold, setThreshold] = useState("");
  const [binColumn, setBinColumn] = useState(null);
  const [bins, setBins] = useState(""); // JSON string for bins
  const [labels, setLabels] = useState(""); // JSON string for labels
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState(null);

  const handleFileChange = async (event) => {
    const selectedFile = event.target.files[0];
    setFile(selectedFile);

    if (selectedFile) {
      const formData = new FormData();
      formData.append("file", selectedFile);

      try {
        const response = await axios.post("http://localhost:8000/preview-columns/", formData, {
          headers: { "Content-Type": "multipart/form-data" },
        });
        setColumns(response.data.columns);
      } catch (error) {
        console.error("Error fetching columns:", error);
      }
    }
  };

  const handleSubmit = async (event) => {
    event.preventDefault();
    setLoading(true);

    const formData = new FormData();
    formData.append("file", file);
    formData.append("is_first_line_header", isFirstLineHeader);
    formData.append("drop_missing", dropMissing);
    if (fillColumn && fillValue) {
      formData.append("fill_missing", JSON.stringify({ [fillColumn]: fillValue }));
    }
    formData.append("categorical_column", categoricalColumn);
    formData.append("threshold", threshold);
    formData.append("bin_column", binColumn);
    formData.append("bins", bins);
    formData.append("labels", labels);

    try {
      const response = await axios.post("http://localhost:8000/upload/", formData, {
        headers: { "Content-Type": "multipart/form-data" },
      });
      setResult(response.data);
    } catch (error) {
      console.error("Error uploading dataset:", error);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div style={{ padding: "2rem" }}>
      <h1>Dataset Processor</h1>
      <form onSubmit={handleSubmit}>
        <div style={{ marginBottom: "1rem" }}>
          <Input
            type="file"
            onChange={handleFileChange}
            fullWidth
            label="Upload Dataset"
          />
        </div>
        <div style={{ marginBottom: "1rem" }}>
          <Checkbox
            isSelected={isFirstLineHeader}
            onChange={(checked) => setIsFirstLineHeader(checked)}
          >
            Columns are in the first line
          </Checkbox>
        </div>
        <div style={{ marginBottom: "1rem" }}>
          <Checkbox
            isSelected={dropMissing}
            onChange={(checked) => setDropMissing(checked)}
          >
            Drop Missing Values
          </Checkbox>
        </div>
        {columns.length > 0 && (
          <div style={{ marginBottom: "1rem" }}>
            <Dropdown>
              <Dropdown.Trigger>
                <Button>{fillColumn || "Select column to fill"}</Button>
              </Dropdown.Trigger>
              <Dropdown.Menu
                selectionMode="single"
                onSelectionChange={(key) => setFillColumn(key)}
              >
                {columns.map((col, index) => (
                  <Dropdown.Item key={col}>{col}</Dropdown.Item>
                ))}
              </Dropdown.Menu>
            </Dropdown>
            {fillColumn && (
              <Input
                label="Fill Value"
                placeholder="Enter value to fill"
                value={fillValue}
                onChange={(e) => setFillValue(e.target.value)}
                fullWidth
              />
            )}
          </div>
        )}
        {columns.length > 0 && (
          <div style={{ marginBottom: "1rem" }}>
            <Dropdown>
              <Dropdown.Trigger>
                <Button>{categoricalColumn || "Select column for categorical"}</Button>
              </Dropdown.Trigger>
              <Dropdown.Menu
                selectionMode="single"
                onSelectionChange={(key) => setCategoricalColumn(key)}
              >
                {columns.map((col, index) => (
                  <Dropdown.Item key={col}>{col}</Dropdown.Item>
                ))}
              </Dropdown.Menu>
            </Dropdown>
          </div>
        )}
        <div style={{ marginBottom: "1rem" }}>
          <Input
            label="Threshold for Categorical Conversion"
            placeholder="Enter threshold"
            value={threshold}
            onChange={(e) => setThreshold(e.target.value)}
            type="number"
            fullWidth
          />
        </div>
        {columns.length > 0 && (
          <div style={{ marginBottom: "1rem" }}>
            <Dropdown>
              <Dropdown.Trigger>
                <Button>{binColumn || "Select column to bin"}</Button>
              </Dropdown.Trigger>
              <Dropdown.Menu
                selectionMode="single"
                onSelectionChange={(key) => setBinColumn(key)}
              >
                {columns.map((col, index) => (
                  <Dropdown.Item key={col}>{col}</Dropdown.Item>
                ))}
              </Dropdown.Menu>
            </Dropdown>
          </div>
        )}
        <div style={{ marginBottom: "1rem" }}>
          <Textarea
            label="Bins (JSON)"
            placeholder="[0, 30, 60, 90]"
            value={bins}
            onChange={(e) => setBins(e.target.value)}
            fullWidth
          />
        </div>
        <div style={{ marginBottom: "1rem" }}>
          <Textarea
            label="Labels for Bins (JSON)"
            placeholder='["Young", "Middle-Aged", "Senior"]'
            value={labels}
            onChange={(e) => setLabels(e.target.value)}
            fullWidth
          />
        </div>
        <div style={{ marginBottom: "1rem" }}>
          <Button type="submit" disabled={!file || loading} auto>
            {loading ? <Loading type="points" /> : "Process Dataset"}
          </Button>
        </div>
      </form>
      {result && (
        <Card style={{ marginTop: "2rem" }}>
          <h3>Processed Dataset</h3>
          <pre>{JSON.stringify(result, null, 2)}</pre>
        </Card>
      )}
    </div>
  );
}
