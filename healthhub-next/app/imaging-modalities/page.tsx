import React from "react";

/**
 * A scientific representation of medical imaging services
 * with structured information and technical styling
 */

// Define service categories with meaningful scientific groupings
const categories = [

  {
    id: "detection",
    name: "Pathology Detection",
    icon: "⚕️",
    description: "Computer vision algorithms specialized in disease and pathology identification"
  },
     {
    id: "diagnostic",
    name: "Diagnostic Analysis",
    icon: "🔍",
    description: "Advanced AI models for diagnostic support and anomaly detection in medical imaging"
  },
  {
    id: "processing",
    name: "Image Processing",
    icon: "⚙️",
    description: "Tools for enhancing, segmenting, and optimizing medical images for analysis"
  },
];

  // Mapping the original services to scientific categories
const serviceCategoryMapping = {
  "Diabetic Retinopathy Detection": "detection",
  "Brain Tumor Segmentation": "diagnostic",
  "Chest Pathology Classification": "detection",
  "COVID-19 Detection": "detection",
  "Image Processing": "processing",
  "Invasive Ductal Carcinoma Detection": "diagnostic",
  "Skin Cancer Detection": "detection",
  "Region of Interest": "processing",
  "Alzheimer's Disease Assessment": "diagnostic",
  // New mappings for added modalities
  "Bone Fracture Analysis": "diagnostic",
  "Liver Lesion Detection": "detection",
  "Cardiac MRI Analysis": "diagnostic",
  "Stroke Lesion Segmentation": "diagnostic",
  "Dental X-Ray Analysis": "detection",
  "Image Denoising": "processing",
  "Prostate Cancer Detection": "detection",
  "Resolution Enhancement": "processing"
};

// Technical capabilities for each service
const technicalCapabilities = {
  "Diabetic Retinopathy Detection": ["Convolutional Neural Networks", "Lesion Segmentation", "Severity Grading"],
  "Brain Tumor Segmentation": ["3D Volumetric Analysis", "Multi-sequence MRI", "Automated Tumor Boundary Detection"],
  "Chest Pathology Classification": ["Multi-label Classification", "Attention Mechanisms", "Localization Maps"],
  "COVID-19 Detection": ["Transfer Learning", "Pneumonia Differentiation", "Severity Assessment"],
  "Image Processing": ["Noise Reduction", "Contrast Enhancement", "Edge Detection"],
  "Invasive Ductal Carcinoma Detection": ["Tissue Classification", "Cell Morphology Analysis", "Histopathological Features"],
  "Skin Cancer Detection": ["Dermoscopic Feature Extraction", "Melanoma vs. Benign Classification", "Lesion Boundary Detection"],
  "Region of Interest": ["Automatic Segmentation", "Multi-scale Analysis", "Feature Extraction"],
  "Alzheimer's Disease Assessment": ["Hippocampal Volume Analysis", "Cortical Thickness Measurement", "Tau/Amyloid Detection"],
  // New capabilities for added modalities
  "Bone Fracture Analysis": ["Fracture Type Classification", "3D Reconstruction", "Displacement Measurement"],
  "Liver Lesion Detection": ["Focal Lesion Characterization", "Texture Analysis", "Volumetric Assessment"],
  "Cardiac MRI Analysis": ["Ventricular Function Assessment", "Myocardial Perfusion", "Wall Motion Analysis"],
  "Stroke Lesion Segmentation": ["Diffusion-Weighted Imaging Analysis", "Penumbra Estimation", "Time Evolution Tracking"],
  "Dental X-Ray Analysis": ["Cavity Detection", "Root Canal Assessment", "Periodontal Disease Classification"],
  "Image Denoising": ["Deep Learning Filters", "Wavelet Decomposition", "Adaptive Thresholding"],
  "Prostate Cancer Detection": ["Multiparametric MRI Analysis", "PI-RADS Classification", "Zonal Segmentation"],
  "Resolution Enhancement": ["Super-Resolution CNNs", "Image Upscaling", "Detail Preservation"]
};

// The main scientific services component
export default function ScientificServices() {
  // Original services data
  const originalServices = [
    {
      title: "Diabetic Retinopathy Detection",
      imagePath: "/modalities_thumbnails/diabetic_retinopathy.jpg",
      description: "Leverage AI to detect diabetic retinopathy with precision.",
      path: "/imaging-modalities/diabetic-retinopathy"
    },
    {
      title: "Brain Tumor Segmentation",
      imagePath: "/modalities_thumbnails/brain_seg.jpg",
      description: "AI-powered 3D MRI segmentation for accurate brain tumor analysis.",
      path: "/imaging-modalities/brain-tumor-segmentation"
    },
    {
      title: "Chest Pathology Classification",
      imagePath: "/modalities_thumbnails/chest.jpg",
      description: "AI-assisted detection of chest pathologies using X-ray imaging.",
      path: "/imaging-modalities/chest-pathology-classification"
    },
    {
      title: "COVID-19 Detection",
      imagePath: "/modalities_thumbnails/covid19.jpeg",
      description: "AI-assisted detection of COVID-19 disease.",
      path: "/imaging-modalities/covid-19-prediction"
    },
    {
      title: "Image Processing",
      imagePath: "/modalities_thumbnails/image-proccesing.jpeg",
      description: "Standard Image Processing Tool.",
      path: "/imaging-modalities/image-processing"
    },
    {
      title: "Invasive Ductal Carcinoma Detection",
      imagePath: "/modalities_thumbnails/breast-cancer.jpg",
      description: "AI-powered detection of invasive ductal carcinoma in breast cancer.",
      path: "/imaging-modalities/invasive-ductal-carcinoma"
    },
    {
      title: "Skin Cancer Detection",
      imagePath: "/modalities_thumbnails/skin-cancer.jpg",
      description: "Detect malignant and benign cancer in dermoscopy images.",
      path: "/imaging-modalities/skin-cancer"
    },
    {
      title: "Region of Interest",
      imagePath: "/modalities_thumbnails/roi.jfif",
      description: "Crop out the region of interest from the image.",
      path: "/imaging-modalities/region-of-interest"
    },
    {
      title: "Alzheimer's Disease Assessment",
      imagePath: "/modalities_thumbnails/alz.png",
      description: "AI-powered detection of Alzheimer's disease biomarkers from brain MRI and PET scans.",
      path: "/imaging-modalities/alzheimers-disease-assessment"
    },
    // New dummy modalities
    {
      title: "Bone Fracture Analysis",
      imagePath: "/modalities_thumbnails/bone.jpg",
      description: "AI-powered detection and classification of bone fractures in X-ray images.",
      path: "/imaging-modalities/bone-fracture-analysis"
    },
    {
      title: "Liver Lesion Detection",
      imagePath: "/modalities_thumbnails/liver-ct.jpg",
      description: "Advanced analysis of CT scans for liver lesion detection and classification.",
      path: "/imaging-modalities/liver-lesion-detection"
    },
    {
      title: "Cardiac MRI Analysis",
      imagePath: "/modalities_thumbnails/cardiac.jpg",
      description: "Comprehensive cardiac function assessment using MRI image analysis.",
      path: "/imaging-modalities/cardiac-mri-analysis"
    },
    {
      title: "Stroke Lesion Segmentation",
      imagePath: "/modalities_thumbnails/stroke.jpg",
      description: "Precise segmentation of ischemic and hemorrhagic stroke lesions in brain MRI.",
      path: "/imaging-modalities/stroke-lesion-segmentation"
    },
    {
      title: "Dental X-Ray Analysis",
      imagePath: "/modalities_thumbnails/dental.jfif",
      description: "AI-assisted diagnosis of dental conditions from panoramic and periapical X-rays.",
      path: "/imaging-modalities/dental-xray-analysis"
    },
    {
      title: "Image Denoising",
      imagePath: "/modalities_thumbnails/denoise.jpg",
      description: "Advanced algorithms for noise reduction in low-quality medical images.",
      path: "/imaging-modalities/image-denoising"
    },
    {
      title: "Prostate Cancer Detection",
      imagePath: "/modalities_thumbnails/healthy_prostate.png",
      description: "AI-based detection and grading of prostate cancer from multiparametric MRI.",
      path: "/imaging-modalities/prostate-cancer-detection"
    },
    {
      title: "Resolution Enhancement",
      imagePath: "/modalities_thumbnails/resolution.png",
      description: "Super-resolution techniques to improve detail in medical imaging.",
      path: "/imaging-modalities/resolution-enhancement"
    }
  ];

  // Enhance services with scientific metadata
  const enhancedServices = originalServices.map(service => ({
    ...service,
    category: serviceCategoryMapping[service.title],
    technicalCapabilities: technicalCapabilities[service.title] || []
  }));

  // Group services by category
  const groupedServices = categories.map(category => ({
    ...category,
    services: enhancedServices.filter(service => service.category === category.id)
  }));

  return (
    <div className="max-w-7xl mx-auto p-6 bg-gray-50">
      {/* Scientific Header */}
      <div className="bg-white p-6 rounded-lg shadow-md mb-8 border-l-4 border-blue-600">
        <div className="flex flex-col md:flex-row justify-between items-start gap-6">
          <div className="flex-1">
            <h1 className="text-3xl font-bold text-gray-800 mb-4">Advanced Medical Imaging Analysis</h1>
            <p className="text-gray-600 mb-4 text-lg">
              Leveraging deep learning architectures and computer vision algorithms to enable
              precision diagnostics and automated medical image analysis.
            </p>
            <div className="flex flex-wrap gap-3 mb-4">
              <span className="bg-blue-100 text-blue-800 px-3 py-1 rounded-full text-sm font-medium">
                Deep Learning Models
              </span>
              <span className="bg-purple-100 text-purple-800 px-3 py-1 rounded-full text-sm font-medium">
                Computer Vision
              </span>
              <span className="bg-green-100 text-green-800 px-3 py-1 rounded-full text-sm font-medium">
                Medical Diagnostics
              </span>
              <span className="bg-red-100 text-red-800 px-3 py-1 rounded-full text-sm font-medium">
                Clinical Validation
              </span>
            </div>
          </div>

          <div className="w-full md:w-64 p-4 bg-blue-50 rounded-lg border border-blue-200">
            <div className="text-center mb-3">
              <span className="inline-block p-2 bg-blue-100 rounded-full">
                <svg className="w-6 h-6 text-blue-700" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 3v2m6-2v2M9 19v2m6-2v2M5 9H3m2 6H3m18-6h-2m2 6h-2M7 19h10a2 2 0 002-2V7a2 2 0 00-2-2H7a2 2 0 00-2 2v10a2 2 0 002 2zM9 9h6v6H9V9z" />
                </svg>
              </span>
            </div>
            <h3 className="text-md font-semibold text-center mb-2">AI Engine</h3>
            <div className="space-y-2 text-sm">
              {/*<div className="flex justify-between">*/}
              {/*  <span className="text-gray-600">CNN Layers</span>*/}
              {/*  <span className="font-medium">18-152</span>*/}
              {/*</div>*/}
              <div className="flex justify-between">
                <span className="text-gray-600">Model Types</span>
                <span className="font-medium">17</span>
              </div>
              <div className="flex justify-between">
                <span className="text-gray-600">Image Formats</span>
                <span className="font-medium">12+</span>
              </div>
              <div className="flex justify-between">
                <span className="text-gray-600">Processing Time</span>
                <span className="font-medium">&lt;500ms</span>
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* Technical Capabilities Overview */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-8">
        {categories.map(category => (
          <div key={category.id} className="bg-white p-4 rounded-lg shadow-md border-l-4 border-blue-500">
            <div className="flex items-center mb-2">
              <span className="text-xl mr-2">{category.icon}</span>
              <h3 className="font-bold text-lg">{category.name}</h3>
            </div>
            <p className="text-sm text-gray-600 mb-2">{category.description}</p>
            <div className="text-xs text-gray-500">
              {category.services?.length || 0} modalities available
            </div>
          </div>
        ))}
      </div>

      {/* Category Sections with Services */}
      {groupedServices.map((category) => (
        <div key={category.id} className="mb-12">
          <div className="flex items-center mb-6">
            <span className="text-2xl mr-2">{category.icon}</span>
            <h2 className="text-2xl font-bold text-gray-800">{category.name}</h2>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-2 gap-6">
            {category.services.map((service) => (
              <div
                key={service.title}
                className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition-shadow border border-gray-200 overflow-hidden"
              >
                <div className="flex flex-col md:flex-row gap-4">
                  <div className="w-full md:w-1/3">
                    <div className="aspect-square overflow-hidden rounded-lg bg-gray-100">
                      <img
                        src={service.imagePath || "/api/placeholder/200/200"}
                        alt={service.title}
                        className="w-full h-full object-cover"
                      />
                    </div>
                  </div>
                  <div className="w-full md:w-2/3">
                    <h3 className="text-lg font-bold text-gray-800 mb-2">{service.title}</h3>
                    <p className="text-sm text-gray-600 mb-3">{service.description}</p>

                    {/* Technical Capabilities */}
                    <div className="mb-4">
                      <p className="text-xs font-medium text-gray-500 mb-1">TECHNICAL CAPABILITIES</p>
                      <div className="flex flex-wrap gap-2">
                        {service.technicalCapabilities.map((capability, index) => (
                          <span key={index} className="inline-block bg-gray-100 text-gray-800 px-2 py-1 rounded text-xs">
                            {capability}
                          </span>
                        ))}
                      </div>
                    </div>

                    {/* Technical Diagram */}
                    <div className="bg-blue-50 p-2 rounded mb-4 hidden md:block">
                      <div className="flex items-center justify-between h-6 px-1">
                        <div className="bg-blue-200 h-2 rounded-full" style={{ width: '60%' }}></div>
                        <span className="text-xs text-blue-700">Input</span>
                      </div>
                      <div className="flex items-center justify-between h-6 px-1">
                        <div className="bg-blue-400 h-2 rounded-full" style={{ width: '75%' }}></div>
                        <span className="text-xs text-blue-700">Processing</span>
                      </div>
                      <div className="flex items-center justify-between h-6 px-1">
                        <div className="bg-blue-600 h-2 rounded-full" style={{ width: '90%' }}></div>
                        <span className="text-xs text-blue-700">Output</span>
                      </div>
                    </div>

                    <a
                      href={service.path}
                      className="inline-block bg-blue-600 hover:bg-blue-700 text-white px-4 py-2 rounded text-sm font-medium transition-colors"
                    >
                      Use service
                    </a>
                  </div>
                </div>
              </div>
            ))}
          </div>
        </div>
      ))}

      {/* Technical Information Footer */}
      <div className="bg-white p-6 rounded-lg shadow-md mt-8 border-t-4 border-gray-800">
        <h3 className="text-lg font-bold mb-3">Technical Information</h3>
        <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
          <div>
            <h4 className="text-md font-semibold mb-2">Model Architecture</h4>
            <ul className="text-sm text-gray-600 space-y-1 list-disc pl-4">
              <li>Custom CNN architectures</li>
              <li>Transfer learning from ImageNet</li>
              <li>Attention mechanisms</li>
              <li>Feature pyramid networks</li>
            </ul>
          </div>
          <div>
            <h4 className="text-md font-semibold mb-2">Technical Standards</h4>
            <ul className="text-sm text-gray-600 space-y-1 list-disc pl-4">
              <li>DICOM compatibility</li>
              <li>HL7 FHIR integration</li>
              <li>HIPAA compliance</li>
              <li>ISO 13485 certification</li>
            </ul>
          </div>
          <div>
            <h4 className="text-md font-semibold mb-2">Deployment Options</h4>
            <ul className="text-sm text-gray-600 space-y-1 list-disc pl-4">
              <li>Cloud-based API</li>
              <li>On-premise deployment</li>
              <li>Edge computing devices</li>
              <li>PACS integration</li>
            </ul>
          </div>
        </div>
        <div className="mt-6 text-sm">
          <p className="text-gray-500">All models undergo rigorous clinical validation and comply with regulatory standards for medical software.</p>
        </div>
      </div>
    </div>
  );
}