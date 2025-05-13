import React from "react";
import { Divider } from "@nextui-org/divider";

/**
 * Icons for scientific visualization
 * Simple SVG icons representing different transcriptomics analysis types
 */
const Icons = {
  CellAnnotation: () => (
    <svg viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg" className="w-full h-full">
      <circle cx="12" cy="12" r="6" stroke="#6366F1" strokeWidth="1.5" fill="none" />
      <circle cx="8" cy="10" r="2" fill="#6366F1" opacity="0.7" />
      <circle cx="14" cy="9" r="1.5" fill="#6366F1" opacity="0.5" />
      <circle cx="15" cy="14" r="2.2" fill="#6366F1" opacity="0.6" />
      <circle cx="10" cy="15" r="1.8" fill="#6366F1" opacity="0.8" />
    </svg>
  ),
  DataConcatenation: () => (
    <svg viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg" className="w-full h-full">
      <rect x="3" y="6" width="7" height="4" rx="1" fill="#6366F1" opacity="0.7" />
      <rect x="14" y="6" width="7" height="4" rx="1" fill="#6366F1" opacity="0.7" />
      <rect x="8.5" y="14" width="7" height="4" rx="1" fill="#6366F1" opacity="0.9" />
      <path d="M10 8H14" stroke="#6366F1" strokeWidth="1.5" strokeLinecap="round" />
      <path d="M12 8V14" stroke="#6366F1" strokeWidth="1.5" strokeLinecap="round" />
    </svg>
  ),
  DataIntegration: () => (
    <svg viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg" className="w-full h-full">
      <rect x="4" y="4" width="6" height="6" rx="1" fill="#6366F1" opacity="0.6" />
      <rect x="4" y="14" width="6" height="6" rx="1" fill="#6366F1" opacity="0.8" />
      <rect x="14" y="9" width="6" height="6" rx="1" fill="#6366F1" opacity="0.7" />
      <path d="M10 7H14V9" stroke="#6366F1" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" />
      <path d="M10 17H12V12H14" stroke="#6366F1" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" />
    </svg>
  ),
  DataPreprocessing: () => (
    <svg viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg" className="w-full h-full">
      <path d="M4 4L20 4" stroke="#6366F1" strokeWidth="1.5" strokeLinecap="round" />
      <path d="M4 8L16 8" stroke="#6366F1" strokeWidth="1.5" strokeLinecap="round" />
      <path d="M4 12L13 12" stroke="#6366F1" strokeWidth="1.5" strokeLinecap="round" />
      <path d="M4 16L10 16" stroke="#6366F1" strokeWidth="1.5" strokeLinecap="round" />
      <path d="M4 20L8 20" stroke="#6366F1" strokeWidth="1.5" strokeLinecap="round" />
      <path d="M20 8L20 20" stroke="#6366F1" strokeWidth="1.5" strokeLinecap="round" strokeDasharray="0.5 2" />
    </svg>
  ),
  ExpressionPlots: () => (
    <svg viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg" className="w-full h-full">
      <path d="M3 20L3 8" stroke="#6366F1" strokeWidth="1.5" strokeLinecap="round" />
      <path d="M3 20L21 20" stroke="#6366F1" strokeWidth="1.5" strokeLinecap="round" />
      <path d="M7 16C7 16 8 8 10 8C12 8 12 14 14 14C16 14 17 4 17 4" stroke="#6366F1" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
      <circle cx="10" cy="8" r="1" fill="#6366F1" />
      <circle cx="14" cy="14" r="1" fill="#6366F1" />
      <circle cx="17" cy="4" r="1" fill="#6366F1" />
    </svg>
  ),
};

/**
 * An array of transcriptomics service objects, each representing a specific analysis type.
 * Each service object contains the following properties:
 *
 * @property {string} title - The title of the analysis type.
 * @property {string} description - A brief description of the analysis type and its capabilities.
 * @property {string} path - The URL path to the detailed page of the analysis.
 * @property {function} icon - React component function that renders an SVG icon
 */
const transcriptomicsServices = [
  {
    title: "Cell Annotation",
    icon: Icons.CellAnnotation,
    description: "AI-powered identification and classification of cell types based on transcriptomic profiles.",
    path: "/transcriptomics/cell-annotation",
  },
  {
    title: "Data Concatenation",
    icon: Icons.DataConcatenation,
    description: "Advanced algorithms for merging multiple transcriptomic datasets while preserving data integrity.",
    path: "/transcriptomics/data-concatenation",
  },
  {
    title: "Data Integration",
    icon: Icons.DataIntegration,
    description: "Multi-modal data integration methods to combine transcriptomic data with other omics platforms.",
    path: "/transcriptomics/data-integration",
  },
  {
    title: "Data Preprocessing",
    icon: Icons.DataPreprocessing,
    description: "State-of-the-art normalization and quality control pipelines for raw transcriptomic data.",
    path: "/transcriptomics/data-preprocessing",
  },
  {
    title: "Expression Plots",
    icon: Icons.ExpressionPlots,
    description: "Interactive visualization tools for gene expression data across cell types and conditions.",
    path: "/transcriptomics/expression-plots",
  },
];

/**
 * ServiceCard component for displaying a service item in the grid
 *
 * @param {Object} props - The props object
 * @param {string} props.title - The title of the service
 * @param {function} props.icon - Icon component function
 * @param {string} props.description - Description of the service
 * @param {string} props.path - URL path to the service details page
 * @returns {JSX.Element} ServiceCard component
 */
function ServiceCard({ title, icon: Icon, description, path }) {
  return (
    <div className="bg-white rounded-lg overflow-hidden shadow-md hover:shadow-lg transition-shadow duration-300 flex flex-col h-full border border-gray-200">
      <div className="h-40 flex items-center justify-center bg-indigo-50 p-6">
        <div className="w-24 h-24">
          <Icon />
        </div>
      </div>
      <div className="p-4 flex-grow flex flex-col">
        <h3 className="text-xl font-semibold text-indigo-900 mb-2">{title}</h3>
        <p className="text-gray-600 text-sm flex-grow">{description}</p>
        <div className="mt-4">
          <a
            href={path}
            className="inline-flex items-center justify-center px-4 py-2 bg-indigo-700 text-white text-sm font-medium rounded hover:bg-indigo-800 transition-colors duration-300 w-full"
          >
            View Details
          </a>
        </div>
      </div>
    </div>
  );
}

/**
 * Renders the Scientific Services component which displays a list of transcriptomics analysis tools.
 *
 * @returns {JSX.Element} The rendered Scientific Services component with integrated ServiceCard.
 */
export default function Services() {
  return (
    <div className="max-w-7xl mx-auto p-6 bg-gray-50">
      <div className="mb-12 pb-6 border-b border-gray-200">
        <h1 className="text-4xl font-bold text-center mb-4 text-indigo-900">
          Transcriptomics Analysis Suite
        </h1>
        <p className="text-lg text-gray-700 text-center max-w-4xl mx-auto">
          Our computational platform provides cutting-edge tools for single-cell and bulk RNA sequencing analysis,
          leveraging machine learning algorithms to extract biological insights from complex transcriptomic datasets.
        </p>
      </div>

      <div className="grid gap-8 grid-cols-1 sm:grid-cols-2 lg:grid-cols-3">
        {transcriptomicsServices.map((service, index) => (
          <ServiceCard
            key={index}
            description={service.description}
            icon={service.icon}
            path={service.path}
            title={service.title}
          />
        ))}
      </div>

      <div className="mt-12 pt-6 border-t border-gray-200">
        <div className="flex items-center justify-center mb-4">
          <span className="h-0.5 w-12 bg-indigo-700 mr-4"></span>
          <h2 className="text-lg font-medium text-indigo-900">Technical Specifications</h2>
          <span className="h-0.5 w-12 bg-indigo-700 ml-4"></span>
        </div>
        <p className="text-sm text-gray-600 text-center max-w-3xl mx-auto">
          All analysis tools implement peer-reviewed computational methods and adhere to established bioinformatics
          standards. Our platform utilizes dimensionality reduction techniques (UMAP, t-SNE, PCA),
          graph-based algorithms, and deep learning models for optimal performance.
        </p>
      </div>
    </div>
  );
}