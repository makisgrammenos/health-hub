import React from "react";
import NextLink from "next/link";

type ServiceCardProps = {
  title: string;
  imagePath: string;
  description: string;
  path: string;
};

/**
 * A card component that displays a service with an image, title, and description.
 * The card is clickable and navigates to the specified path when clicked.
 *
 * @component
 * @param {Object} props - The properties object.
 * @param {string} props.title - The title of the service.
 * @param {string} props.imagePath - The path to the image representing the service.
 * @param {string} props.description - A brief description of the service.
 * @param {string} props.path - The path to navigate to when the card is clicked.
 * @returns {JSX.Element} The rendered ServiceCard component.
 */
const ServiceCard: React.FC<ServiceCardProps> = ({
  title,
  imagePath,
  description,
  path,
}) => {
  return (
    <div className="group relative cursor-pointer rounded-lg overflow-hidden shadow-lg hover:shadow-2xl transition-shadow">
      <NextLink href={path}>
        <div className="relative w-full h-60">
          <img
            src={imagePath}
            alt={title}
            className="object-cover w-full h-full transition-transform duration-500 ease-in-out group-hover:scale-105 group-hover:blur-sm"
          />
        </div>
        <div className="absolute inset-0 bg-black bg-opacity-50 group-hover:bg-opacity-30 transition-all"></div>
        <div className="absolute bottom-4 left-4">
          <h2 className="text-lg font-bold text-white">{title}</h2>
          <p className="text-sm text-gray-300">{description}</p>
        </div>
      </NextLink>
    </div>
  );
};

export default ServiceCard;
