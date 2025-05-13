'use client';

import React, { useState, useRef, useEffect } from "react";

export default function ROISelector() {
  const [uploadedImage, setUploadedImage] = useState(null);
  const [base64Image, setBase64Image] = useState(null);
  const [drawing, setDrawing] = useState(false);
  const [roi, setRoi] = useState({});
  const [croppedBase64, setCroppedBase64] = useState(null);
  const [croppedFilename, setCroppedFilename] = useState(null); // New state for cropped filename
  const [uploadedFilename, setUploadedFilename] = useState(null); // New state for uploaded filename
  const [isEditing, setIsEditing] = useState(false); // New state for editing
  const [editMode, setEditMode] = useState(null); // 'move' or 'resize'
  const [currentHandle, setCurrentHandle] = useState(null); // e.g., 'top-left'
  const canvasRef = useRef(null);
  const imageRef = useRef(null);
  const websocketRef = useRef(null);
  const roiRef = useRef({}); // Ref to store current ROI and editing info

  // Maximum canvas dimensions
  const MAX_WIDTH = 800;
  const MAX_HEIGHT = 600;

  const HANDLE_SIZE = 8;

  const getHandles = (roi) => {
    const { x, y, width, height } = roi;
    return [
      { position: 'top-left', x: x, y: y },
      { position: 'top-center', x: x + width / 2, y: y },
      { position: 'top-right', x: x + width, y: y },
      { position: 'middle-left', x: x, y: y + height / 2 },
      { position: 'middle-right', x: x + width, y: y + height / 2 },
      { position: 'bottom-left', x: x, y: y + height },
      { position: 'bottom-center', x: x + width / 2, y: y + height },
      { position: 'bottom-right', x: x + width, y: y + height },
    ];
  };

  const drawCanvas = () => {
    if (!canvasRef.current || !imageRef.current) return;
    const ctx = canvasRef.current.getContext("2d");
    ctx.clearRect(0, 0, canvasRef.current.width, canvasRef.current.height);
    ctx.drawImage(imageRef.current, 0, 0, canvasRef.current.width, canvasRef.current.height);

    if (roi.width && roi.height) {
      ctx.strokeStyle = "rgba(255, 0, 0, 0.7)";
      ctx.lineWidth = 3;
      ctx.setLineDash([6]);
      ctx.strokeRect(roi.x, roi.y, roi.width, roi.height);
      ctx.setLineDash([]);

      // Draw handles
      const handles = getHandles(roi);
      handles.forEach(handle => {
        ctx.fillStyle = "rgba(255, 0, 0, 0.7)";
        ctx.fillRect(handle.x - HANDLE_SIZE / 2, handle.y - HANDLE_SIZE / 2, HANDLE_SIZE, HANDLE_SIZE);
      });
    }
  };

  // Establish WebSocket connection once when component mounts
  useEffect(() => {
    websocketRef.current = new WebSocket("ws://localhost:8000/imaging/roi/ws");

    websocketRef.current.onopen = () => {
      console.log("WebSocket connection established");
    };

    websocketRef.current.onmessage = (event) => {
      const data = JSON.parse(event.data);
      console.log("Message from server:", data);

      if (data.message === "Image uploaded successfully" && data.filename) {
        setUploadedFilename(data.filename); // Store the uploaded filename
      }

      if (data.message === "Image cropped successfully" && data.image_base64 && data.cropped_filename) {
        // Display the cropped image next to the canvas in real time
        setCroppedBase64(data.image_base64);
        setCroppedFilename(data.cropped_filename); // Store the cropped filename
      }

      if (data.error) {
        console.error("Server Error:", data.error);
      }
    };

    websocketRef.current.onerror = (error) => {
      console.error("WebSocket error:", error);
    };

    websocketRef.current.onclose = (event) => {
      console.log("WebSocket connection closed:", event.reason);
      // Optionally implement reconnection logic here
    };

    // Cleanup WebSocket connection when component unmounts
    return () => {
      if (websocketRef.current) {
        websocketRef.current.close();
      }
    };
  }, []);

  // Redraw canvas whenever ROI or uploadedImage changes
  useEffect(() => {
    drawCanvas();
  }, [roi, uploadedImage]);

  const handleUpload = async (event) => {
    const file = event.target.files[0];
    if (!file) return;
    const imageUrl = URL.createObjectURL(file);
    setUploadedImage(imageUrl);

    // Convert image to base64
    const reader = new FileReader();
    reader.onload = () => {
      const result = reader.result;
      // Strip off "data:image/...;base64," prefix
      const base64Data = result.split(",")[1];
      setBase64Image(base64Data);
      // Once we have the base64 image, send it via WebSocket
      if (websocketRef.current && websocketRef.current.readyState === WebSocket.OPEN) {
        websocketRef.current.send(JSON.stringify({
          action: "upload_image",
          image: base64Data,
        }));
      } else {
        console.error("WebSocket is not open. Unable to send image.");
      }
    };
    reader.readAsDataURL(file);

    const img = new Image();
    img.src = imageUrl;
    img.onload = () => {
      imageRef.current = img;

      const scale = Math.min(MAX_WIDTH / img.width, MAX_HEIGHT / img.height, 1);
      const scaledWidth = img.width * scale;
      const scaledHeight = img.height * scale;

      const canvas = canvasRef.current;
      canvas.width = scaledWidth;
      canvas.height = scaledHeight;

      const ctx = canvas.getContext("2d");
      ctx.drawImage(img, 0, 0, scaledWidth, scaledHeight);
    };
  };

  const handleMouseDown = (e) => {
    if (!canvasRef.current) return;
    const rect = canvasRef.current.getBoundingClientRect();
    const startX = e.clientX - rect.left;
    const startY = e.clientY - rect.top;

    if (roi.width && roi.height) {
      const handles = getHandles(roi);
      // Check if the mouse is over any handle
      const handle = handles.find(h => 
        Math.abs(h.x - startX) <= HANDLE_SIZE / 2 && Math.abs(h.y - startY) <= HANDLE_SIZE / 2
      );

      if (handle) {
        // Start resizing
        setIsEditing(true);
        setEditMode('resize');
        setCurrentHandle(handle.position);
        roiRef.current = {
          original: { ...roi },
          startX,
          startY,
        };
        return;
      }

      // Check if the mouse is inside the ROI for moving
      if (
        startX >= roi.x &&
        startX <= roi.x + roi.width &&
        startY >= roi.y &&
        startY <= roi.y + roi.height
      ) {
        // Start moving
        setIsEditing(true);
        setEditMode('move');
        roiRef.current = {
          original: { ...roi },
          startX,
          startY,
        };
        return;
      }
    }

    // If not editing, start drawing a new ROI
    setDrawing(true);
    const newRoi = {
      x: startX,
      y: startY,
      width: 0,
      height: 0,
    };
    setRoi(newRoi);
    roiRef.current = newRoi;
  };

  const handleMouseMove = (e) => {
    if (!canvasRef.current || !imageRef.current) return;
    const rect = canvasRef.current.getBoundingClientRect();
    const currentX = e.clientX - rect.left;
    const currentY = e.clientY - rect.top;

    if (drawing) {
      let newX = roiRef.current.x;
      let newY = roiRef.current.y;
      let newWidth = currentX - roiRef.current.x;
      let newHeight = currentY - roiRef.current.y;

      // Normalize the ROI to ensure width and height are positive
      if (newWidth < 0) {
        newX = currentX;
        newWidth = Math.abs(newWidth);
      }
      if (newHeight < 0) {
        newY = currentY;
        newHeight = Math.abs(newHeight);
      }

      const updatedRoi = {
        x: newX,
        y: newY,
        width: newWidth,
        height: newHeight,
      };

      setRoi(updatedRoi);
      roiRef.current = updatedRoi;

      const ctx = canvasRef.current.getContext("2d");
      ctx.clearRect(0, 0, canvasRef.current.width, canvasRef.current.height);
      ctx.drawImage(imageRef.current, 0, 0, canvasRef.current.width, canvasRef.current.height);

      ctx.strokeStyle = "rgba(255, 0, 0, 0.7)";
      ctx.lineWidth = 3;
      ctx.setLineDash([6]);
      ctx.strokeRect(updatedRoi.x, updatedRoi.y, updatedRoi.width, updatedRoi.height);
      ctx.setLineDash([]);

      // Draw handles
      const handles = getHandles(updatedRoi);
      handles.forEach(handle => {
        ctx.fillStyle = "rgba(255, 0, 0, 0.7)";
        ctx.fillRect(handle.x - HANDLE_SIZE / 2, handle.y - HANDLE_SIZE / 2, HANDLE_SIZE, HANDLE_SIZE);
      });
    } else if (isEditing) {
      const ctx = canvasRef.current.getContext("2d");
      ctx.clearRect(0, 0, canvasRef.current.width, canvasRef.current.height);
      ctx.drawImage(imageRef.current, 0, 0, canvasRef.current.width, canvasRef.current.height);

      if (editMode === 'move') {
        const deltaX = currentX - roiRef.current.startX;
        const deltaY = currentY - roiRef.current.startY;
        const newX = roiRef.current.original.x + deltaX;
        const newY = roiRef.current.original.y + deltaY;

        const updatedRoi = {
          ...roiRef.current.original,
          x: newX,
          y: newY,
        };
        setRoi(updatedRoi);
        roiRef.current = updatedRoi;
      } else if (editMode === 'resize' && currentHandle) {
        let { x, y, width, height } = roiRef.current.original;
        switch (currentHandle) {
          case 'top-left':
            width += x - currentX;
            height += y - currentY;
            x = currentX;
            y = currentY;
            break;
          case 'top-center':
            height += y - currentY;
            y = currentY;
            break;
          case 'top-right':
            width = currentX - x;
            height += y - currentY;
            y = currentY;
            break;
          case 'middle-left':
            width += x - currentX;
            x = currentX;
            break;
          case 'middle-right':
            width = currentX - x;
            break;
          case 'bottom-left':
            width += x - currentX;
            height = currentY - y;
            x = currentX;
            break;
          case 'bottom-center':
            height = currentY - y;
            break;
          case 'bottom-right':
            width = currentX - x;
            height = currentY - y;
            break;
          default:
            break;
        }

        // Prevent negative width and height
        if (width < 0) {
          x += width;
          width = Math.abs(width);
        }
        if (height < 0) {
          y += height;
          height = Math.abs(height);
        }

        const updatedRoi = { x, y, width, height };
        setRoi(updatedRoi);
        roiRef.current = updatedRoi;
      }

      // Redraw the ROI and handles
      ctx.strokeStyle = "rgba(255, 0, 0, 0.7)";
      ctx.lineWidth = 3;
      ctx.setLineDash([6]);
      ctx.strokeRect(roiRef.current.x, roiRef.current.y, roiRef.current.width, roiRef.current.height);
      ctx.setLineDash([]);

      // Draw handles
      const handles = getHandles(roiRef.current);
      handles.forEach(handle => {
        ctx.fillStyle = "rgba(255, 0, 0, 0.7)";
        ctx.fillRect(handle.x - HANDLE_SIZE / 2, handle.y - HANDLE_SIZE / 2, HANDLE_SIZE, HANDLE_SIZE);
      });
    }
  };

  const handleMouseUp = () => {
    if (drawing) {
      setDrawing(false);

      // After drawing, send the crop action if the websocket is open
      if (
        websocketRef.current &&
        websocketRef.current.readyState === WebSocket.OPEN &&
        roiRef.current.width > 0 &&
        roiRef.current.height > 0
      ) {
        const scale = Math.min(
          MAX_WIDTH / imageRef.current.width,
          MAX_HEIGHT / imageRef.current.height,
          1
        );

        const unscaledRoi = {
          x: Math.round(roiRef.current.x / scale),
          y: Math.round(roiRef.current.y / scale),
          width: Math.round(roiRef.current.width / scale),
          height: Math.round(roiRef.current.height / scale),
        };

        websocketRef.current.send(
          JSON.stringify({
            action: "crop",
            roi: unscaledRoi,
            filename: uploadedFilename, // Send the filename of the uploaded image
          })
        );
      }
    } else if (isEditing) {
      setIsEditing(false);
      setEditMode(null);
      setCurrentHandle(null);

      // After editing, send the updated ROI to the backend
      if (
        websocketRef.current &&
        websocketRef.current.readyState === WebSocket.OPEN &&
        roiRef.current.width > 0 &&
        roiRef.current.height > 0
      ) {
        const scale = Math.min(
          MAX_WIDTH / imageRef.current.width,
          MAX_HEIGHT / imageRef.current.height,
          1
        );

        const unscaledRoi = {
          x: Math.round(roiRef.current.x / scale),
          y: Math.round(roiRef.current.y / scale),
          width: Math.round(roiRef.current.width / scale),
          height: Math.round(roiRef.current.height / scale),
        };

        websocketRef.current.send(
          JSON.stringify({
            action: "crop",
            roi: unscaledRoi,
            filename: uploadedFilename, // Send the filename of the uploaded image
          })
        );
      }
    }
  };

  const handleDelete = () => {
    setRoi({});
    roiRef.current = {};
    setCroppedBase64(null);
    setCroppedFilename(null);
    const ctx = canvasRef.current.getContext("2d");
    ctx.clearRect(0, 0, canvasRef.current.width, canvasRef.current.height);
    ctx.drawImage(imageRef.current, 0, 0, canvasRef.current.width, canvasRef.current.height);
  };

  const handleDownload = () => {
    if (!croppedFilename) {
      alert("No cropped image available for download.");
      return;
    }

    // Construct the download URL
    const downloadUrl = `http://localhost:8000/imaging/roi/download/${croppedFilename}`;

    // Create a temporary link to trigger the download
    const link = document.createElement("a");
    link.href = downloadUrl;
    link.download = croppedFilename; // Set the desired file name
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };

  return (
    <div className="container mx-auto p-4">
      <h1 className="text-4xl font-bold mb-6 text-center">Region of Interest Selector</h1>

      <div className="flex flex-col items-center mb-6">
        <input
          type="file"
          accept="image/*"
          onChange={handleUpload}
          className="mb-4 text-sm text-gray-700 file:py-2 file:px-4 file:rounded file:bg-blue-50 file:text-blue-700 file:border file:border-blue-500 hover:file:bg-blue-100"
        />
      </div>

      {uploadedImage && (
        <div className="flex flex-col md:flex-row items-start space-x-6">
          {/* Image and Canvas */}
          <div className="relative">
            <canvas
              ref={canvasRef}
              style={{
                border: "2px solid #ccc",
                borderRadius: "8px",
                background: "white",
                maxWidth: "100%",
                boxShadow: "0px 4px 10px rgba(0, 0, 0, 0.1)",
                cursor: "crosshair",
              }}
              onMouseDown={handleMouseDown}
              onMouseMove={handleMouseMove}
              onMouseUp={handleMouseUp}
            />
          </div>

          {/* ROI Info Panel (Show while drawing or when a ROI exists) */}
          {(drawing || (roi.width && roi.height)) && (
            <div className="bg-gray-100 p-4 rounded shadow-sm text-sm min-w-[200px] mt-4 md:mt-0">
              <h2 className="font-bold mb-2">Current ROI</h2>
              <p>X: {Math.round(roi.x) || 0}</p>
              <p>Y: {Math.round(roi.y) || 0}</p>
              <p>Width: {Math.round(roi.width) || 0}</p>
              <p>Height: {Math.round(roi.height) || 0}</p>
            </div>
          )}

          {/* Cropped Image Preview and Download/Delete Buttons */}
          {croppedBase64 && (
            <div className="bg-white p-4 rounded shadow-sm text-sm max-w-[300px] mt-4 md:mt-0">
              <h2 className="font-bold mb-2">Selected Region of Interest Preview</h2>
              <img
                src={`data:image/jpeg;base64,${croppedBase64}`}
                alt="Cropped Preview"
                style={{ maxWidth: "100%" }}
              />
              <div className="flex space-x-2 mt-4">
                <button
                  onClick={handleDownload}
                  className="px-4 py-2 bg-blue-500 text-white rounded hover:bg-blue-600"
                >
                  Download Cropped ROI
                </button>
                {/*<button*/}
                {/*  onClick={handleDownload}*/}
                {/*  className="px-4 py-2 bg-blue-500 text-white rounded hover:bg-blue-600"*/}
                {/*>*/}
                {/*  Show statistics*/}
                {/*</button>*/}
                <button
                  onClick={handleDelete}
                  className="px-4 py-2 bg-red-500 text-white rounded hover:bg-red-600"
                >
                  Delete ROI
                </button>
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
