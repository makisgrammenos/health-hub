'use client'
import { Button } from "@nextui-org/button";
import { Image } from "@nextui-org/image";
import { Card, CardBody, CardFooter } from "@nextui-org/card";
import { Divider, Chip, Badge } from "@nextui-org/react";
import {
  Microscope,
  Activity,
  Database,
  Award,
  Users,
  ChevronRight,
  FileText,
  Search,
  AlertCircle,
  Brain,
  Dna,
  ImageIcon,
  HeartPulse,
  FlaskConical
} from "lucide-react";

type ToolCardProps = {
  title: string;
  description: string;
  buttonText: string;
  buttonLink: string;
  imageSrc: string;
  icon: React.ReactNode;
  toolCount?: number;
  isNew?: boolean;
  tags?: string[];
};

// Feature stat item component
const StatItem = ({ icon, value, label }: { icon: React.ReactNode, value: string, label: string }) => (
  <div className="flex flex-col items-center p-4 bg-white/80 backdrop-blur-sm rounded-lg shadow-sm">
    <div className="mb-2 text-blue-600">{icon}</div>
    <p className="text-2xl font-bold text-gray-800">{value}</p>
    <p className="text-xs text-gray-500 mt-1">{label}</p>
  </div>
);

// Partner logo component
const PartnerLogo = ({ src, alt }: { src: string, alt: string }) => (
  <div className="bg-white p-3 rounded-lg shadow-sm hover:shadow-md transition-shadow flex items-center justify-center h-16">
    <Image
      src={src}
      alt={alt}
      width={120}
      height={40}
      className="object-contain max-h-10"
    />
  </div>
);

const ToolCard: React.FC<ToolCardProps> = ({
  title,
  description,
  buttonText,
  buttonLink,
  imageSrc,
  icon,
  toolCount,
  isNew = false,
  tags = [],
}) => {
  return (
    <Card className="border border-transparent hover:border-blue-200 transition-all duration-300 overflow-hidden">
      <div className="relative">
        <Image
          src={imageSrc}
          alt={title}
          width={400}
          height={250}
          className="object-cover w-full h-48"
        />
        <div className="absolute inset-0 bg-gradient-to-t from-black/70 via-transparent to-transparent" />

        {isNew && (
          <Badge
            content="NEW"
            color="primary"
            placement="top-right"
            className="absolute top-2 right-2 z-10"
          />
        )}

        <div className="absolute bottom-0 left-0 p-4 text-white">
          <div className="flex items-center gap-2 mb-1">
            {icon}
            <h3 className="text-lg font-bold">{title}</h3>
          </div>
          {toolCount && (
            <div className="flex items-center gap-1 text-xs">
              <Database size={12} />
              <span>{toolCount} tools available</span>
            </div>
          )}
        </div>
      </div>

      <CardBody className="p-4">
        <p className="text-sm text-gray-600">{description}</p>

        {tags.length > 0 && (
          <div className="flex flex-wrap gap-1 mt-3">
            {tags.map((tag, index) => (
              <Chip key={index} size="sm" variant="flat" color="primary" className="text-xs">
                {tag}
              </Chip>
            ))}
          </div>
        )}
      </CardBody>

      <CardFooter className="pt-0">
        <Button
          as="a"
          href={buttonLink}
          color="primary"
          variant="flat"
          radius="lg"
          endContent={<ChevronRight size={16} />}
          className="w-full"
        >
          {buttonText}
        </Button>
      </CardFooter>
    </Card>
  );
};

// Scientific publication card component
const PublicationCard = ({
  title,
  journal,
  date,
  impact
}: {
  title: string,
  journal: string,
  date: string,
  impact: number
}) => (
  <Card className="bg-white shadow-sm hover:shadow-md transition-all">
    <CardBody className="p-4">
      <div className="flex justify-between items-start mb-2">
        <Chip size="sm" color="primary" variant="flat">{journal}</Chip>
        <Chip size="sm" variant="flat">IF: {impact}</Chip>
      </div>
      <h4 className="text-sm font-semibold line-clamp-2 mb-2">{title}</h4>
      <div className="flex justify-between items-center text-xs text-gray-500">
        <span>{date}</span>
        <Button size="sm" variant="light" color="primary" className="p-1">
          View paper
        </Button>
      </div>
    </CardBody>
  </Card>
);

export default function Home() {
  // Mock data for statistics
  const stats = [
    { value: "14+", label: "PARTNERED INSTITUTIONS", icon: <Users size={24} /> },
    { value: "12", label: "AI TOOLS DEVELOPED", icon: <Brain size={24} /> },
    { value: "9.2M+", label: "PATIENT RECORDS ANALYZED", icon: <Database size={24} /> },
    { value: "19", label: "CASE STUDIES", icon: <FileText size={24} /> }
  ];

  // Mock data for latest publications
  const publications = [
    {
      title: "Deep Learning Models for Early Detection of Diabetic Retinopathy: A Multi-center Study",
      journal: "Nature Medicine",
      date: "May 2025",
      impact: 53.44
    },
    {
      title: "AI-driven Analysis of Genomic Data for Personalized Cancer Treatment",
      journal: "Cell",
      date: "April 2025",
      impact: 41.58
    },
    {
      title: "Integration of Multi-modal Health Data: A European Perspective",
      journal: "The Lancet Digital Health",
      date: "March 2025",
      impact: 21.83
    }
  ];

  return (
    <main className="w-full min-h-screen bg-gray-50 text-gray-900">
      {/* Hero Section with Video Background */}
      <section className="relative w-full h-[70vh] flex flex-col items-center justify-center text-center px-4 overflow-hidden">
        {/* Video Background */}
        <video
          autoPlay
          loop
          muted
          playsInline
          className="absolute inset-0 w-full h-full object-cover z-0"
        >
          <source src="hero/video.webm" type="video/webm" />
          Your browser does not support the video tag.
        </video>

        {/* Gradient Overlay */}
        <div className="absolute inset-0 bg-gradient-to-b from-black/70 via-black/60 to-black/30 z-10" />

        {/* Hero Content */}
        <div className="relative z-20 max-w-5xl mx-auto">
          <div className="flex items-center justify-center mb-4">
            <Chip color="primary" variant="dot" className="text-white"> Digital Innovation Hub</Chip>
          </div>

          <h1 className="text-4xl md:text-6xl font-extrabold tracking-tight text-white mb-6">
            HealthHub<span className="text-blue-400"> AI </span>Platform
          </h1>

          <p className="mt-4 text-lg md:text-xl font-light max-w-3xl mx-auto text-gray-100">
            Accelerating the Digital Transformation of Healthcare & Pharmaceuticals through
            <span className="text-blue-300 font-medium"> Advanced AI Applications</span>
          </p>

          <div className="flex flex-wrap gap-4 justify-center mt-8">
            <Button
              href="#tools"
              as="a"
              size="lg"
              color="primary"
              variant="shadow"
              radius="full"
              className="bg-gradient-to-r from-blue-600 to-indigo-600 font-medium px-8"
            >
              Explore AI Tools
            </Button>

            <Button
              href="#about"
              as="a"
              size="lg"
              variant="bordered"
              radius="full"
              className="text-white border-white/60 font-medium hover:bg-white/10"
            >
              Learn More
            </Button>
          </div>
        </div>
      </section>

      {/* Key Statistics Section */}
      <section className="py-12 px-6">
        <div className="max-w-7xl mx-auto">
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
            {stats.map((stat, index) => (
              <StatItem key={index} icon={stat.icon} value={stat.value} label={stat.label} />
            ))}
          </div>
        </div>
      </section>

      {/* About Section with Gradient Background */}
      <section
        id="about"
        className="py-20 bg-gradient-to-br from-blue-900 via-indigo-900 to-blue-800 text-white"
      >
        <div className="max-w-7xl mx-auto px-6">
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-16 items-center">
            <div>
              <Chip color="primary" variant="flat" className="mb-4">About HealthHub</Chip>
              <h2 className="text-3xl md:text-4xl font-bold mb-6">
                Transforming European Healthcare Through Advanced AI Solutions
              </h2>
              <p className="text-lg text-blue-100 mb-6">
                HealthHub is a European Digital Innovation Hub dedicated to accelerating the adoption of
                artificial intelligence in healthcare and pharmaceutical sectors across Europe.
              </p>

              <div className="space-y-4 mb-8">
                <div className="flex items-start gap-3">
                  <div className="mt-1 text-blue-300">
                    <Award size={20} />
                  </div>
                  <div>
                    <h3 className="font-semibold text-lg">Excellence in Research</h3>
                    <p className="text-blue-100 text-sm">
                      Our consortium brings together leading research institutions, hospitals, and industrial partners
                      to develop and validate cutting-edge AI technologies.
                    </p>
                  </div>
                </div>

                <div className="flex items-start gap-3">
                  <div className="mt-1 text-blue-300">
                    <Activity size={20} />
                  </div>
                  {/*<div>*/}
                  {/*  <h3 className="font-semibold text-lg">Clinical Validation</h3>*/}
                  {/*  <p className="text-blue-100 text-sm">*/}
                  {/*    All tools undergo rigorous clinical validation in real-world healthcare settings across multiple*/}
                  {/*    European centers to ensure reliability and effectiveness.*/}
                  {/*  </p>*/}
                  {/*</div>*/}
                </div>

                <div className="flex items-start gap-3">
                  <div className="mt-1 text-blue-300">
                    <AlertCircle size={20} />
                  </div>
                  <div>
                    <h3 className="font-semibold text-lg">Ethical AI Development</h3>
                    <p className="text-blue-100 text-sm">
                      We adhere to the highest ethical standards in AI development, ensuring privacy, transparency,
                      fairness, and responsible use of health data.
                    </p>
                  </div>
                </div>
              </div>

              <Button
                as="a"
                href="/about"
                color="primary"
                variant="shadow"
                radius="lg"
                className="bg-white text-blue-800 hover:bg-blue-50"
              >
                Learn About Our Consortium
              </Button>
            </div>

            <div className="relative hidden lg:block">
              <div className="absolute -top-16 -right-16 w-64 h-64 bg-blue-400/10 rounded-full blur-3xl" />
              <div className="absolute -bottom-8 -left-8 w-48 h-48 bg-indigo-400/10 rounded-full blur-2xl" />

              <div className="relative bg-gradient-to-br from-white/10 to-white/5 backdrop-blur-sm p-8 rounded-2xl border border-white/10">
                <div className="mb-8">
                  <h3 className="text-xl font-semibold mb-4">Key Focus Areas</h3>
                  <div className="grid grid-cols-2 gap-3">
                    <div className="flex items-center gap-2 bg-white/10 p-3 rounded-lg">
                      <ImageIcon size={20} className="text-blue-300" />
                      <span className="text-sm">Medical Imaging</span>
                    </div>
                    <div className="flex items-center gap-2 bg-white/10 p-3 rounded-lg">
                      <Dna size={20} className="text-blue-300" />
                      <span className="text-sm">Genomics</span>
                    </div>
                    <div className="flex items-center gap-2 bg-white/10 p-3 rounded-lg">
                      <HeartPulse size={20} className="text-blue-300" />
                      <span className="text-sm">Clinical Informatics</span>
                    </div>
                    <div className="flex items-center gap-2 bg-white/10 p-3 rounded-lg">
                      <Brain size={20} className="text-blue-300" />
                      <span className="text-sm">Neuroscience</span>
                    </div>
                  </div>
                </div>

                {/*<div>*/}
                {/*  <h3 className="text-xl font-semibold mb-4">Funding Sources</h3>*/}
                {/*  <div className="flex flex-wrap gap-3">*/}
                {/*    <Chip size="sm" variant="flat">Horizon Europe</Chip>*/}
                {/*    <Chip size="sm" variant="flat">Digital Europe Programme</Chip>*/}
                {/*    <Chip size="sm" variant="flat">European Research Council</Chip>*/}
                {/*  </div>*/}
                {/*</div>*/}
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* Tools Section with Enhanced Cards */}
      <section
        id="tools"
        className="py-20 px-6 md:px-12 bg-gray-50"
      >
        <div className="max-w-7xl mx-auto">
          <div className="text-center mb-12">
            <Chip color="primary" variant="flat" className="mb-2">AI-Powered Solutions</Chip>
            <h2 className="text-3xl md:text-4xl font-bold text-gray-800 mb-4">
              Advanced AI Tools for Healthcare
            </h2>
            <p className="text-lg text-gray-600 max-w-3xl mx-auto">
              Access our suite of clinically-validated AI tools designed to enhance diagnostic
              accuracy, optimize treatment pathways, and accelerate research.
            </p>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8">
            <ToolCard
              title="Medical Imaging"
              description="AI-powered tools for advanced image analysis across multiple modalities including MRI, CT, X-ray, and fundus photography."
              buttonText="Explore Imaging Tools"
              buttonLink="/imaging-modalities"
              imageSrc="/banner/medical-imaging.jpg"
              icon={<ImageIcon size={20} />}
              toolCount={12}
              isNew={true}
              tags={["Diabetic Retinopathy", "Chest X-Ray", "Mammography", "Brain MRI"]}
            />

            <ToolCard
              title="Clinical Informatics"
              description="Tools for EHR data analysis, clinical decision support, risk stratification, and patient outcome prediction."
              buttonText="Explore Clinical Tools"
              buttonLink="/clinical-tools"
              imageSrc="/banner/clinical.jpg"
              icon={<Activity size={20} />}
              toolCount={8}
              tags={["Risk Prediction", "Disease Progression", "Treatment Response"]}
            />

            <ToolCard
              title="Multi-omics Analysis"
              description="Advanced tools for genomics, proteomics, metabolomics, and integrated multi-omics data analysis and interpretation."
              buttonText="Explore Molecular Tools"
              buttonLink="/transcriptomics"
              imageSrc="/banner/molecular.png"
              icon={<Dna size={20} />}
              toolCount={7}
              tags={["Transcriptomics", "Precision Medicine", "Single-cell Analysis"]}
            />
          </div>

          <div className="mt-12 text-center">
            <Button
              as="a"
              href="/all-tools"
              color="primary"
              variant="light"
              endContent={<Search size={16} />}
              size="lg"
            >
              View All Available Tools
            </Button>
          </div>
        </div>
      </section>

      {/*/!* Featured Research Section *!/*/}
      {/*<section className="py-20 px-6 bg-gradient-to-br from-gray-50 to-blue-50">*/}
      {/*  <div className="max-w-7xl mx-auto">*/}
      {/*    <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-12">*/}
      {/*      <div>*/}
      {/*        <Chip color="primary" variant="flat" className="mb-2">Latest Research</Chip>*/}
      {/*        <h2 className="text-3xl font-bold text-gray-800">Featured Publications</h2>*/}
      {/*      </div>*/}
      {/*      <Button*/}
      {/*        as="a"*/}
      {/*        href="/publications"*/}
      {/*        color="primary"*/}
      {/*        variant="light"*/}
      {/*        className="mt-4 md:mt-0"*/}
      {/*      >*/}
      {/*        View All Publications*/}
      {/*      </Button>*/}
      {/*    </div>*/}

      {/*    <div className="grid grid-cols-1 md:grid-cols-3 gap-6">*/}
      {/*      {publications.map((pub, index) => (*/}
      {/*        <PublicationCard*/}
      {/*          key={index}*/}
      {/*          title={pub.title}*/}
      {/*          journal={pub.journal}*/}
      {/*          date={pub.date}*/}
      {/*          impact={pub.impact}*/}
      {/*        />*/}
      {/*      ))}*/}
      {/*    </div>*/}
      {/*  </div>*/}
      {/*</section>*/}

      {/* Partners Section */}
      <section className="py-20 px-6 bg-white">
        <div className="max-w-7xl mx-auto">
          <div className="text-center mb-12">
            <h2 className="text-2xl font-bold text-gray-800 mb-2">Our Partners & Collaborators</h2>
            <p className="text-gray-600">Working together to advance healthcare through AI innovation</p>
          </div>

          <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-6 gap-4">
            {[...Array(6)].map((_, index) => (
              <PartnerLogo
                key={index}
                src={`/api/placeholder/120/60`}
                alt={`Partner ${index + 1}`}
              />
            ))}
          </div>
        </div>
      </section>

      {/* CTA Section */}
      <section className="py-16 px-6 bg-gradient-to-r from-blue-700 to-indigo-800 text-white">
        <div className="max-w-5xl mx-auto text-center">
          <h2 className="text-3xl font-bold mb-4">Ready to Transform Your Research?</h2>
          <p className="text-lg text-blue-100 mb-8 max-w-3xl mx-auto">
            Join leading medical institutions across Europe using our AI tools to accelerate discovery
            and improve patient outcomes.
          </p>
          <div className="flex flex-wrap gap-4 justify-center">
            <Button
              as="a"
              href="/register"
              color="primary"
              size="lg"
              radius="full"
              className="bg-white text-blue-700 font-medium px-8"
            >
              Register Your Institution
            </Button>
            <Button
              as="a"
              href="/contact"
              variant="bordered"
              size="lg"
              radius="full"
              className="text-white border-white"
            >
              Contact Our Team
            </Button>
          </div>
        </div>
      </section>

      {/* Footer */}
      <footer className="bg-gray-900 text-gray-400 py-12 px-6">
        <div className="max-w-7xl mx-auto">
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-8">
            <div>
              <h3 className="text-white text-lg font-bold mb-4">HealthHub</h3>
              <p className="text-sm mb-4">
                European Digital Innovation Hub for the Digital Transformation of Health & Pharmaceuticals
                through Artificial Intelligence Applications
              </p>
              <div>
                <Chip size="sm" variant="flat" color="primary">EU Funded Project</Chip>
              </div>
            </div>

            <div>
              <h4 className="text-white text-base font-semibold mb-4">Resources</h4>
              <ul className="space-y-2 text-sm">
                <li><a href="/tools" className="hover:text-white transition-colors">AI Tools</a></li>
                <li><a href="/publications" className="hover:text-white transition-colors">Publications</a></li>
                <li><a href="/case-studies" className="hover:text-white transition-colors">Case Studies</a></li>
                <li><a href="/research" className="hover:text-white transition-colors">Research Areas</a></li>
              </ul>
            </div>

            <div>
              <h4 className="text-white text-base font-semibold mb-4">About</h4>
              <ul className="space-y-2 text-sm">
                <li><a href="/consortium" className="hover:text-white transition-colors">Consortium</a></li>
                <li><a href="/partners" className="hover:text-white transition-colors">Partners</a></li>
                <li><a href="/funding" className="hover:text-white transition-colors">Funding</a></li>
                <li><a href="/contact" className="hover:text-white transition-colors">Contact Us</a></li>
              </ul>
            </div>

            <div>
              <h4 className="text-white text-base font-semibold mb-4">Legal</h4>
              <ul className="space-y-2 text-sm">
                <li><a href="/privacy" className="hover:text-white transition-colors">Privacy Policy</a></li>
                <li><a href="/terms" className="hover:text-white transition-colors">Terms of Use</a></li>
                <li><a href="/data-policy" className="hover:text-white transition-colors">Data Policy</a></li>
                <li><a href="/accessibility" className="hover:text-white transition-colors">Accessibility</a></li>
              </ul>
            </div>
          </div>

          <Divider className="my-8 bg-gray-800" />

          <div className="flex flex-col md:flex-row justify-between items-center">
            <p className="text-sm">
              © {new Date().getFullYear()} HealthHub Consortium. All rights reserved.
            </p>
            <p className="text-xs mt-2 md:mt-0">
              Funded by the European Union under grant agreement No. 123456.
            </p>
          </div>
        </div>
      </footer>
    </main>
  );
}