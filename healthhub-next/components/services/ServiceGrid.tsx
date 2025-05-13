import React from "react";
import ServiceCard from "../../components/services/ServiceCard";
import { Divider } from "@nextui-org/divider";

interface Service {
  title: string;
  imagePath: string;
  description: string;
  path: string;
}

interface ServiceGridProps {
  services: Service[];
  title: string;
  description?: string;
}

const ServiceGrid: React.FC<ServiceGridProps> = ({ services, title, description }) => {
  return (
    <div className="max-w-7xl mx-auto p-6">
      <h1 className="text-4xl font-extrabold text-center mb-8">{title}</h1>
      {description && (
        <p className="text-lg text-gray-600 text-center mb-12">{description}</p>
      )}
      <div className="grid gap-6 grid-cols-1 sm:grid-cols-2 md:grid-cols-3">
        {services.map((service, index) => (
          <ServiceCard
            key={index}
            description={service.description}
            imagePath={service.imagePath}
            path={service.path}
            title={service.title}
          />
        ))}
      </div>
      <Divider className="my-4" />
    </div>
  );
};

export default ServiceGrid;
