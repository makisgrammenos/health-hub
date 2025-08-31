'use client'
import React, { useState } from "react";

/**
 * Medical Imaging Analysis Platform
 * Professional clinical interface for medical professionals and researchers
 */

const ScientificServices = () => {
  const [selectedModality, setSelectedModality] = useState("all");
  const [viewMode, setViewMode] = useState("grid"); // grid or table

  // Clinical modality categories
  const modalities = [
    { id: "all", name: "All Modalities", code: "ALL" },
    { id: "radiology", name: "Radiology", code: "RAD" },
    { id: "mri", name: "MRI", code: "MRI" },
    { id: "pathology", name: "Pathology", code: "PATH" },
    { id: "ophthalmology", name: "Ophthalmology", code: "OPTH" },
    { id: "nuclear", name: "Nuclear Medicine", code: "NM" },
    { id: "dermatology", name: "Dermatology", code: "DERM" },
    { id: "processing", name: "Processing", code: "PROC" }
  ];

  // Clinical service data
  const services = [
    {
      id: "DR-001",
      title: "Diabetic Retinopathy Detection",
      modality: "ophthalmology",
      imagePath: "/modalities_thumbnails/diabetic_retinopathy.jpg",
      description: "Automated grading of diabetic retinopathy severity from fundus photographs using deep learning classification",
      path: "/imaging-modalities/diabetic-retinopathy",
      specifications: {
        sensitivity: "95.2%",
        specificity: "93.8%",
        auc: "0.974",
        processingTime: "280ms",
        validation: "10,542 images"
      },
      technical: {
        algorithm: "ResNet-50",
        input: "Color fundus photography (45° FOV)",
        output: "ICDR severity scale (0-4)",
        standards: "DICOM, HL7 FHIR"
      }
    },
    {
      id: "BT-001",
      title: "Brain Tumor Segmentation",
      modality: "mri",
      imagePath: "/modalities_thumbnails/brain_seg.jpg",
      description: "3D volumetric segmentation of gliomas in multimodal MRI with tumor subregion delineation",
      path: "/imaging-modalities/brain-tumor-segmentation",
      specifications: {
        sensitivity: "94.8%",
        specificity: "95.1%",
        dice: "0.892",
        processingTime: "420ms",
        validation: "BraTS 2021 dataset"
      },
      technical: {
        algorithm: "3D U-Net",
        input: "T1, T1c, T2, FLAIR sequences",
        output: "Segmentation masks (WT, TC, ET)",
        standards: "DICOM, NIfTI"
      }
    },
    {
      id: "CX-001",
      title: "Chest Pathology Classification",
      modality: "radiology",
      imagePath: "/modalities_thumbnails/chest.jpg",
      description: "Multi-label classification system for 14 thoracic pathologies from frontal chest radiographs",
      path: "/imaging-modalities/chest-pathology-classification",
      specifications: {
        sensitivity: "92.7%",
        specificity: "94.2%",
        auc: "0.941",
        processingTime: "190ms",
        validation: "ChestX-ray14 dataset"
      },
      technical: {
        algorithm: "DenseNet-121",
        input: "Chest X-ray (PA/AP view)",
        output: "Probability scores for 14 pathologies",
        standards: "DICOM"
      }
    },
    {
      id: "CV-001",
      title: "COVID-19 Detection",
      modality: "radiology",
      imagePath: "/modalities_thumbnails/covid19.jpeg",
      description: "CT-based COVID-19 pneumonia detection with severity scoring and lung involvement quantification",
      path: "/imaging-modalities/covid-19-prediction",
      specifications: {
        sensitivity: "96.1%",
        specificity: "95.3%",
        auc: "0.982",
        processingTime: "350ms",
        validation: "COVIDx CT-2 dataset"
      },
      technical: {
        algorithm: "COVIDNet-CT",
        input: "Chest CT (HRCT)",
        output: "COVID probability + CO-RADS score",
        standards: "DICOM"
      }
    },
    {
      id: "IP-001",
      title: "Image Processing Suite",
      modality: "processing",
      imagePath: "/modalities_thumbnails/image-proccesing.jpeg",
      description: "Medical image preprocessing toolkit for enhancement, normalization, and artifact reduction",
      path: "/imaging-modalities/image-processing",
      specifications: {
        snr: "+12.4 dB",
        specificity: "N/A",
        sensitivity: "N/A",
        processingTime: "150ms",
        validation: "Multi-modal"
      },
      technical: {
        algorithm: "Adaptive filters",
        input: "Any medical image format",
        output: "Enhanced image",
        standards: "DICOM, NIfTI, NRRD"
      }
    },
    {
      id: "BC-001",
      title: "Breast Density Estimation",
      modality: "Mammography",
      imagePath: "/modalities_thumbnails/mammograph.jpg",
      description: "Automated breast density estimation from mammographic images based on BI-RADS A–D categories.",
      path: "/imaging-modalities/breast-density",
      specifications: {
        sensitivity: "94.1%",
        specificity: "93.7%",
        auc: "0.96",
        processingTime: "120ms/image",
        validation: "Mayo Clinic dataset"
      },
      technical: {
        algorithm: "Inception-v3",
        input: "Mammography images (JPEG, PNG, DICOM)",
        output: "Breast density classification (BI-RADS A, B, C, D)",
        standards: "DICOM"
      }
    },
    {
      id: "SK-001",
      title: "Skin Cancer Detection",
      modality: "dermatology",
      imagePath: "/modalities_thumbnails/skin-cancer.jpg",
      description: "Melanoma vs benign lesion classification using dermoscopic image analysis",
      path: "/imaging-modalities/skin-cancer",
      specifications: {
        sensitivity: "91.8%",
        specificity: "89.2%",
        auc: "0.937",
        processingTime: "220ms",
        validation: "ISIC 2020 dataset"
      },
      technical: {
        algorithm: "EfficientNet-B4",
        input: "Dermoscopy images",
        output: "Malignancy probability",
        standards: "DICOM"
      }
    },
    {
      id: "ROI-001",
      title: "Region of Interest Extraction",
      modality: "processing",
      imagePath: "/modalities_thumbnails/roi.jfif",
      description: "Automated anatomical structure localization and extraction for targeted analysis",
      path: "/imaging-modalities/region-of-interest",
      specifications: {
        sensitivity: "94.3%",
        specificity: "92.7%",
        iou: "0.89",
        processingTime: "80ms",
        validation: "Multi-organ"
      },
      technical: {
        algorithm: "YOLOv5",
        input: "Any medical image",
        output: "Bounding boxes + masks",
        standards: "DICOM"
      }
    },
    {
      id: "AD-001",
      title: "Alzheimer's Disease Assessment",
      modality: "nuclear",
      imagePath: "/modalities_thumbnails/alz.png",
      description: "Quantitative assessment of amyloid and tau burden from PET imaging with MCI progression risk",
      path: "/imaging-modalities/alzheimers-disease-assessment",
      specifications: {
        sensitivity: "89.3%",
        specificity: "87.6%",
        auc: "0.912",
        processingTime: "520ms",
        validation: "ADNI cohort"
      },
      technical: {
        algorithm: "3D CNN + SVM",
        input: "FDG-PET, Amyloid-PET",
        output: "SUVr + risk score",
        standards: "DICOM, NIfTI"
      }
    },
    {
      id: "FX-001",
      title: "Bone Fracture Analysis",
      modality: "radiology",
      imagePath: "/modalities_thumbnails/bone.jpg",
      description: "Automated detection and AO/OTA classification of fractures in musculoskeletal radiographs",
      path: "/imaging-modalities/bone-fracture-analysis",
      specifications: {
        sensitivity: "94.2%",
        specificity: "93.1%",
        auc: "0.961",
        processingTime: "240ms",
        validation: "MURA dataset"
      },
      technical: {
        algorithm: "RetinaNet",
        input: "X-ray (multi-view)",
        output: "Fracture location + type",
        standards: "DICOM"
      }
    },
    {
      id: "LV-001",
      title: "Liver Lesion Detection",
      modality: "radiology",
      imagePath: "/modalities_thumbnails/liver-ct.jpg",
      description: "Multi-phase CT analysis for hepatocellular carcinoma detection with LI-RADS classification",
      path: "/imaging-modalities/liver-lesion-detection",
      specifications: {
        sensitivity: "92.9%",
        specificity: "91.4%",
        auc: "0.948",
        processingTime: "410ms",
        validation: "LiTS dataset"
      },
      technical: {
        algorithm: "nn-UNet",
        input: "Multi-phase CT",
        output: "LI-RADS score (1-5)",
        standards: "DICOM"
      }
    },
    {
      id: "CR-001",
      title: "Cardiac MRI Analysis",
      modality: "mri",
      imagePath: "/modalities_thumbnails/cardiac.jpg",
      description: "Automated cardiac function assessment including ejection fraction and wall motion scoring",
      path: "/imaging-modalities/cardiac-mri-analysis",
      specifications: {
        sensitivity: "95.7%",
        specificity: "94.1%",
        correlation: "r=0.97",
        processingTime: "480ms",
        validation: "UK Biobank"
      },
      technical: {
        algorithm: "4D CNN-LSTM",
        input: "Cine SSFP sequences",
        output: "EF%, EDV, ESV, strain",
        standards: "DICOM"
      }
    },
    {
      id: "ST-001",
      title: "Stroke Lesion Segmentation",
      modality: "mri",
      imagePath: "/modalities_thumbnails/stroke.jpg",
      description: "Acute ischemic stroke core and penumbra segmentation from diffusion and perfusion MRI",
      path: "/imaging-modalities/stroke-lesion-segmentation",
      specifications: {
        sensitivity: "93.1%",
        specificity: "94.7%",
        dice: "0.871",
        processingTime: "390ms",
        validation: "ISLES 2022"
      },
      technical: {
        algorithm: "DeepMedic",
        input: "DWI, PWI, ADC maps",
        output: "Core + penumbra volumes",
        standards: "DICOM, NIfTI"
      }
    },
    {
      id: "DN-001",
      title: "Dental X-Ray Analysis",
      modality: "radiology",
      imagePath: "/modalities_thumbnails/dental.jfif",
      description: "Comprehensive dental pathology detection including caries, periodontal disease, and periapical lesions",
      path: "/imaging-modalities/dental-xray-analysis",
      specifications: {
        sensitivity: "90.5%",
        specificity: "88.9%",
        auc: "0.923",
        processingTime: "180ms",
        validation: "Clinical dataset"
      },
      technical: {
        algorithm: "Mask R-CNN",
        input: "Panoramic/Periapical X-ray",
        output: "Pathology annotations",
        standards: "DICOM"
      }
    },
    {
      id: "DN-002",
      title: "Medical Image Denoising",
      modality: "processing",
      imagePath: "/modalities_thumbnails/denoise.jpg",
      description: "Deep learning-based noise reduction for low-dose imaging protocols",
      path: "/imaging-modalities/image-denoising",
      specifications: {
        sensitivity: "N/A",
        specificity: "N/A",
        psnr: "+8.2 dB",
        processingTime: "120ms",
        validation: "Multi-modal"
      },
      technical: {
        algorithm: "DnCNN",
        input: "Noisy medical image",
        output: "Denoised image",
        standards: "DICOM"
      }
    },
    {
      id: "PC-001",
      title: "Prostate Cancer Detection",
      modality: "mri",
      imagePath: "/modalities_thumbnails/healthy_prostate.png",
      description: "Multiparametric MRI analysis for prostate cancer detection with PI-RADS v2.1 scoring",
      path: "/imaging-modalities/prostate-cancer-detection",
      specifications: {
        sensitivity: "91.4%",
        specificity: "89.7%",
        auc: "0.934",
        processingTime: "460ms",
        validation: "PROSTATEx"
      },
      technical: {
        algorithm: "3D ResNet",
        input: "T2, DWI, DCE-MRI",
        output: "PI-RADS score + map",
        standards: "DICOM"
      }
    },
    {
      id: "SR-001",
      title: "Resolution Enhancement",
      modality: "processing",
      imagePath: "/modalities_thumbnails/resolution.png",
      description: "Super-resolution reconstruction for improved spatial resolution in medical imaging",
      path: "/imaging-modalities/resolution-enhancement",
      specifications: {
        sensitivity: "N/A",
        specificity: "N/A",
        upscaling: "4x",
        processingTime: "200ms",
        validation: "Multi-modal"
      },
      technical: {
        algorithm: "ESRGAN",
        input: "Low-resolution image",
        output: "SR image (4x)",
        standards: "DICOM"
      }
    }
  ];

  // Filter services
  const filteredServices = selectedModality === "all" 
    ? services 
    : services.filter(s => s.modality === selectedModality);

  // Count by modality
  const getModalityCount = (modalityId) => {
    return modalityId === "all" ? services.length : services.filter(s => s.modality === modalityId).length;
  };

  return (
    <div className="min-h-screen bg-gray-50">
      {/* Clinical Header */}
      <header className="bg-white border-b-2 border-gray-300">
        <div className="max-w-7xl mx-auto px-4 py-4">
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-4">
              <div className="w-10 h-10 bg-blue-900 rounded flex items-center justify-center">
                <svg className="w-6 h-6 text-white" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 3v2m6-2v2M9 19v2m6-2v2M5 9H3m2 6H3m18-6h-2m2 6h-2M7 19h10a2 2 0 002-2V7a2 2 0 00-2-2H7a2 2 0 00-2 2v10a2 2 0 002 2zM9 9h6v6H9V9z" />
                </svg>
              </div>
              <div>
                <h1 className="text-2xl font-bold text-gray-900">
                  Medical Imaging Analysis Platform
                </h1>
                <p className="text-sm text-gray-600">Clinical Decision Support System v4.2.1</p>
              </div>
            </div>
            
            <div className="flex items-center gap-4">
              <div className="text-right">
                <div className="text-xs text-gray-500 uppercase tracking-wide">System Status</div>
                <div className="flex items-center gap-2">
                  <span className="w-2 h-2 bg-green-500 rounded-full"></span>
                  <span className="text-sm font-medium text-gray-700">Operational</span>
                </div>
              </div>
              <div className="text-right">
                <div className="text-xs text-gray-500 uppercase tracking-wide">Compliance</div>
                <div className="text-sm font-medium text-gray-700">HIPAA • ISO 13485</div>
              </div>
            </div>
          </div>
        </div>
      </header>

      {/* Key Metrics Bar */}
      <div className="bg-blue-900 text-white">
        <div className="max-w-7xl mx-auto px-4 py-3">
          <div className="grid grid-cols-6 gap-4 text-center">
            <div>
              <div className="text-xl font-bold">2,421,847</div>
              <div className="text-xs opacity-80">Total Analyses</div>
            </div>
            <div>
              <div className="text-xl font-bold">93.5%</div>
              <div className="text-xs opacity-80">Mean Accuracy</div>
            </div>
            <div>
              <div className="text-xl font-bold">147</div>
              <div className="text-xs opacity-80">Institutions</div>
            </div>
            <div>
              <div className="text-xl font-bold">17</div>
              <div className="text-xs opacity-80">AI Models</div>
            </div>
            <div>
              <div className="text-xl font-bold">&lt;500ms</div>
              <div className="text-xs opacity-80">Avg. Processing</div>
            </div>
            <div>
              <div className="text-xl font-bold">99.97%</div>
              <div className="text-xs opacity-80">Uptime SLA</div>
            </div>
          </div>
        </div>
      </div>

      {/* Control Panel */}
      <div className="bg-white border-b border-gray-200">
        <div className="max-w-7xl mx-auto px-4 py-4">
          <div className="flex items-center justify-between">
            {/* Modality Filter */}
            <div className="flex items-center gap-3">
              <span className="text-sm font-medium text-gray-700">Filter by Modality:</span>
              <div className="flex gap-2">
                {modalities.map(mod => (
                  <button
                    key={mod.id}
                    onClick={() => setSelectedModality(mod.id)}
                    className={`px-3 py-1.5 text-sm font-medium border rounded transition-colors ${
                      selectedModality === mod.id
                        ? 'bg-blue-900 text-white border-blue-900'
                        : 'bg-white text-gray-700 border-gray-300 hover:bg-gray-50'
                    }`}
                  >
                    {mod.code} ({getModalityCount(mod.id)})
                  </button>
                ))}
              </div>
            </div>

            {/* View Toggle */}
            <div className="flex items-center gap-2">
              <span className="text-sm font-medium text-gray-700">View:</span>
              <button
                onClick={() => setViewMode('grid')}
                className={`px-3 py-1.5 text-sm font-medium border rounded-l transition-colors ${
                  viewMode === 'grid'
                    ? 'bg-gray-800 text-white border-gray-800'
                    : 'bg-white text-gray-700 border-gray-300 hover:bg-gray-50'
                }`}
              >
                Grid
              </button>
              <button
                onClick={() => setViewMode('table')}
                className={`px-3 py-1.5 text-sm font-medium border rounded-r transition-colors ${
                  viewMode === 'table'
                    ? 'bg-gray-800 text-white border-gray-800'
                    : 'bg-white text-gray-700 border-gray-300 hover:bg-gray-50'
                }`}
              >
                Table
              </button>
            </div>
          </div>
        </div>
      </div>

      {/* Main Content Area */}
      <div className="max-w-7xl mx-auto px-4 py-6">
        {viewMode === 'grid' ? (
          // Grid View
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {filteredServices.map(service => (
              <div key={service.id} className="bg-white border border-gray-300 rounded-lg overflow-hidden hover:shadow-lg transition-shadow">
                {/* Image */}
                <div className="h-40 bg-gray-100 relative">
                  <img
                    src={service.imagePath || "/api/placeholder/400/200"}
                    alt={service.title}
                    className="w-full h-full object-cover"
                    onError={(e) => {
                      e.target.style.display = 'none';
                      e.target.parentElement.classList.add('flex', 'items-center', 'justify-center');
                      e.target.parentElement.innerHTML = `
                        <div class="text-center">
                          <div class="text-3xl text-gray-400 mb-2">📊</div>
                          <div class="text-xs text-gray-500 uppercase tracking-wide">No Preview</div>
                        </div>
                      `;
                    }}
                  />
                  <div className="absolute top-2 left-2 bg-black bg-opacity-75 text-white px-2 py-1 text-xs font-mono">
                    {service.id}
                  </div>
                </div>

                {/* Content */}
                <div className="p-4">
                  <h3 className="font-bold text-gray-900 mb-2">{service.title}</h3>
                  <p className="text-sm text-gray-600 mb-4 leading-relaxed">
                    {service.description}
                  </p>

                  {/* Specifications */}
                  <div className="border-t border-gray-200 pt-3 mb-3">
                    <div className="grid grid-cols-2 gap-2 text-xs">
                      <div>
                        <span className="text-gray-500">Sensitivity:</span>
                        <span className="ml-1 font-semibold text-gray-900">{service.specifications.sensitivity}</span>
                      </div>
                      <div>
                        <span className="text-gray-500">Specificity:</span>
                        <span className="ml-1 font-semibold text-gray-900">{service.specifications.specificity}</span>
                      </div>
                      <div>
                        <span className="text-gray-500">Processing:</span>
                        <span className="ml-1 font-semibold text-gray-900">{service.specifications.processingTime}</span>
                      </div>
                      <div>
                        <span className="text-gray-500">Validation:</span>
                        <span className="ml-1 font-semibold text-gray-900">{service.specifications.validation}</span>
                      </div>
                    </div>
                  </div>

                  {/* Technical Info */}
                  <div className="bg-gray-50 -mx-4 -mb-4 px-4 py-3 border-t border-gray-200">
                    <div className="text-xs space-y-1">
                      <div>
                        <span className="text-gray-500">Algorithm:</span>
                        <span className="ml-1 font-mono text-gray-700">{service.technical.algorithm}</span>
                      </div>
                      <div>
                        <span className="text-gray-500">Standards:</span>
                        <span className="ml-1 text-gray-700">{service.technical.standards}</span>
                      </div>
                    </div>
                    
                    <a
                      href={service.path}
                      className="mt-3 block w-full bg-blue-900 hover:bg-blue-800 text-white text-center py-2 rounded text-sm font-medium transition-colors"
                    >
                      Open Analysis Module →
                    </a>
                  </div>
                </div>
              </div>
            ))}
          </div>
        ) : (
          // Table View
          <div className="bg-white border border-gray-300 rounded-lg overflow-hidden">
            <table className="w-full">
              <thead className="bg-gray-100 border-b border-gray-300">
                <tr>
                  <th className="px-4 py-3 text-left text-xs font-medium text-gray-700 uppercase tracking-wider">ID</th>
                  <th className="px-4 py-3 text-left text-xs font-medium text-gray-700 uppercase tracking-wider">Service</th>
                  <th className="px-4 py-3 text-left text-xs font-medium text-gray-700 uppercase tracking-wider">Modality</th>
                  <th className="px-4 py-3 text-left text-xs font-medium text-gray-700 uppercase tracking-wider">Performance</th>
                  <th className="px-4 py-3 text-left text-xs font-medium text-gray-700 uppercase tracking-wider">Algorithm</th>
                  <th className="px-4 py-3 text-left text-xs font-medium text-gray-700 uppercase tracking-wider">Processing</th>
                  <th className="px-4 py-3 text-left text-xs font-medium text-gray-700 uppercase tracking-wider">Action</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-200">
                {filteredServices.map(service => (
                  <tr key={service.id} className="hover:bg-gray-50">
                    <td className="px-4 py-3 text-sm font-mono text-gray-900">{service.id}</td>
                    <td className="px-4 py-3">
                      <div>
                        <div className="text-sm font-medium text-gray-900">{service.title}</div>
                        <div className="text-xs text-gray-500">{service.technical.input}</div>
                      </div>
                    </td>
                    <td className="px-4 py-3 text-sm text-gray-700">
                      {modalities.find(m => m.id === service.modality)?.code}
                    </td>
                    <td className="px-4 py-3">
                      <div className="text-xs">
                        <div>Sens: {service.specifications.sensitivity}</div>
                        <div>Spec: {service.specifications.specificity}</div>
                      </div>
                    </td>
                    <td className="px-4 py-3 text-sm font-mono text-gray-700">
                      {service.technical.algorithm}
                    </td>
                    <td className="px-4 py-3 text-sm text-gray-700">
                      {service.specifications.processingTime}
                    </td>
                    <td className="px-4 py-3">
                      <a
                        href={service.path}
                        className="text-blue-900 hover:text-blue-700 text-sm font-medium"
                      >
                        Launch →
                      </a>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        )}
      </div>

      {/* Technical Documentation Footer */}
      <footer className="bg-gray-800 text-gray-300 mt-12">
        <div className="max-w-7xl mx-auto px-4 py-8">
          <div className="grid grid-cols-4 gap-8">
            <div>
              <h4 className="text-white font-semibold mb-3">Technical Specifications</h4>
              <ul className="text-sm space-y-1">
                <li>• NVIDIA DGX A100 Infrastructure</li>
                <li>• TensorRT Optimization</li>
                <li>• REST API v2.0</li>
                <li>• WebSocket Streaming</li>
              </ul>
            </div>
            <div>
              <h4 className="text-white font-semibold mb-3">Regulatory Compliance</h4>
              <ul className="text-sm space-y-1">
                <li>• FDA 510(k) Pending</li>
                <li>• HIPAA Compliant</li>
                <li>• ISO 13485:2016</li>
                <li>• CE Mark (MDR)</li>
              </ul>
            </div>
            <div>
              <h4 className="text-white font-semibold mb-3">Integration Standards</h4>
              <ul className="text-sm space-y-1">
                <li>• DICOM 3.0</li>
                <li>• HL7 FHIR R4</li>
                <li>• IHE Profiles</li>
                <li>• PACS/RIS Compatible</li>
              </ul>
            </div>
            <div>
              <h4 className="text-white font-semibold mb-3">Clinical Validation</h4>
              <ul className="text-sm space-y-1">
                <li>• 47 Peer-reviewed Publications</li>
                <li>• 12 Clinical Trials</li>
                <li>• 2.4M+ Images Analyzed</li>
                <li>• Multi-center Validation</li>
              </ul>
            </div>
          </div>
          
          <div className="border-t border-gray-700 mt-6 pt-6 text-center text-sm">
            <p>© 2025 Medical Imaging Analysis Platform • Research and Clinical Decision Support System</p>
            <p className="mt-2 text-xs text-gray-500">
              This system is intended for use by qualified healthcare professionals. Not for diagnostic use without clinical correlation.
            </p>
          </div>
        </div>
      </footer>
    </div>
  );
};

export default ScientificServices;