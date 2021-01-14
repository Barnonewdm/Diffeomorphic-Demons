/**
 * This is a small tool that shows how to use the diffeomorphic demons algorithm.
 * The user can choose if diffeomorphic, additive or compositive demons should be used.
 * The user can also choose the type of demons forces, or other parameters;
 *
 * \author Tom Vercauteren, INRIA & Mauna Kea Technologies
 */


#include <itkCommand.h>
#include <itkDiffeomorphicDemonsRegistrationFilter.h>
#include <itkDisplacementFieldJacobianDeterminantFilter.h>
#include <itkFastSymmetricForcesDemonsRegistrationFilter.h>
#include <itkGridForwardWarpImageFilter.h>
#include <itkHistogramMatchingImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkMultiResolutionPDEDeformableRegistration.h>
#include <itkTransformFileReader.h>
#include <itkTransformToDeformationFieldSource.h>
#include <itkVectorCentralDifferenceImageFunction.h>
#include <itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h>
#include <itkWarpHarmonicEnergyCalculator.h>
#include <itkWarpImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include <metaCommand.h>

#include <errno.h>
#include <iostream>
#include <limits.h>


struct arguments
{
  // std::string  fixedImageFile;  /* -f option */
  std::string  movingImageFile; /* -m option */
  std::string  inputFieldFile;  /* -b option */
  // std::string  inputTransformFile;  /* -p option */
  std::string  outputImageFile; /* -o option */
  // std::string  outputFieldFile; /* -O option */
  // std::string  trueFieldFile;   /* -r option */
  // std::vector<unsigned int> numIterations;   /* -i option */
  // float sigmaDef;               /* -s option */
  // float sigmaUp;                /* -g option */
  // float maxStepLength;          /* -l option */
  // unsigned int updateRule;      /* -a option */
  // unsigned int gradientType;    /* -t option */
  // bool useHistogramMatching;    /* -e option */
  // unsigned int verbosity;       /* -d option */
	bool useNearestInterpolation; /* -I option */

  friend std::ostream& operator<< (std::ostream& o, const arguments& args)
    {
		return o;
		}
};

void help_callback()
{
  std::cout<<std::endl;
	// std::cout<<"DemonsWarp -m input.hdr -o output.hdr -b field.mha [-I] (no interpolation)"<<std::endl;
  exit( EXIT_FAILURE );
};

int atoi_check( const char * str )
{
  return 0;
}


std::vector<unsigned int> parseUIntVector( const std::string & str)
{
  std::vector<unsigned int> vect;

  std::string::size_type crosspos = str.find('x',0);

  if (crosspos == std::string::npos)
    {
    // only one uint
    vect.push_back( static_cast<unsigned int>( atoi_check(str.c_str()) ));
    return vect;
    }

  // first uint
  vect.push_back( static_cast<unsigned int>(
                     atoi_check( (str.substr(0,crosspos)).c_str()  ) ));

  while(true)
    {
    std::string::size_type crossposfrom = crosspos;
    crosspos =  str.find('x',crossposfrom+1);

    if (crosspos == std::string::npos)
      {
      vect.push_back( static_cast<unsigned int>(
                         atoi_check( (str.substr(crossposfrom+1,str.length()-crossposfrom-1)).c_str()  ) ));
      return vect;
      }

    vect.push_back( static_cast<unsigned int>(
                       atoi_check( (str.substr(crossposfrom+1,crosspos-crossposfrom-1)).c_str()  ) ));
    }
}


void parseOpts (int argc, char **argv, struct arguments & args)
{
  // Command line parser
  MetaCommand command;
  command.SetParseFailureOnUnrecognizedOption( true );
  command.SetHelpCallBack(help_callback);

  // Fill some information about the software
  command.SetAuthor("Tom Vercauteren");

  command.SetAcknowledgments("This work stems from the author's CIFRE PhD thesis at INRIA (Asclepios team) and Mauna Kea Technologies");

  command.SetDescription("Basic image registration tool with the diffeomorphic demons algorithm.");

  // Define parsing options
  /* command.SetOption("FixedImageFile","f",true,"Fixed image filename");
  command.SetOptionLongTag("FixedImageFile","fixed-image");
  command.AddOptionField("FixedImageFile","filename",MetaCommand::STRING,true); */

  command.SetOption("MovingImageFile","m",true,"Moving image filename");
  command.SetOptionLongTag("MovingImageFile","moving-image");
  command.AddOptionField("MovingImageFile","filename",MetaCommand::STRING,true);

  command.SetOption("InputFieldFile","b",false,"Input field filename");
  command.SetOptionLongTag("InputFieldFile","input-field");
  command.AddOptionField("InputFieldFile","filename",MetaCommand::STRING,true);

  /* command.SetOption("InputTransformFile","p",false,"Input transform filename");
  command.SetOptionLongTag("InputTransformFile","input-transform");
  command.AddOptionField("InputTransformFile","filename",MetaCommand::STRING,true); */

  command.SetOption("OutputImageFile","o",false,"Output image filename");
  command.SetOptionLongTag("OutputImageFile","output-image");
  command.AddOptionField("OutputImageFile","filename",MetaCommand::STRING,true,"output.mha");

	command.SetOption("UseNearestInterpolation","I",false,"Use NO interpolation");
  command.SetOptionLongTag("UseNearestInterpolation","use-nearest-interpolation");
  command.AddOptionField("UseNearestInterpolation","boolval",MetaCommand::FLAG,false);

  /* command.SetOption("OutputFieldFile","O",false,"Generate the output field and optionally specify a filename");
  command.SetOptionLongTag("OutputFieldFile","output-field");
  command.AddOptionField("OutputFieldFile","filename",MetaCommand::STRING,false,"OUTPUTIMAGENAME-field.mha");

  command.SetOption("TrueFieldFile","r",false,"Specify a \"true\" field to compare the registration result with (useful for synthetic experiments)");
  command.SetOptionLongTag("TrueFieldFile","true-field");
  command.AddOptionField("TrueFieldFile","filename",MetaCommand::STRING,true);

  command.SetOption("NumberOfIterationsPerLevels","i",false,"List of number of iterations for each multi-scale pyramid level < UINTx...xUINT >");
  command.SetOptionLongTag("NumberOfIterationsPerLevels","num-iterations");
  command.AddOptionField("NumberOfIterationsPerLevels","uintvect",MetaCommand::STRING,true,"15x10x5");

  command.SetOption("DeformationFieldSigma","s",false,"Smoothing sigma for the deformation field (pixel units). Setting it value below 0.5 means no smoothing will be performed");
  command.SetOptionLongTag("DeformationFieldSigma","def-field-sigma");
  command.AddOptionField("DeformationFieldSigma","floatval",MetaCommand::FLOAT,true,"1.5");

  command.SetOption("UpdateFieldSigma","g",false,"Smoothing sigma for the update field (pixel units). Setting it below 0.5 means no smoothing will be performed");
  command.SetOptionLongTag("UpdateFieldSigma","up-field-sigma");
  command.AddOptionField("UpdateFieldSigma","floatval",MetaCommand::FLOAT,true,"0.0");

  command.SetOption("MaximumUpdateStepLength","l",false,"Maximum length of an update vector (pixel units). Setting it to 0 implies no restrictions will be made on the step length");
  command.SetOptionLongTag("MaximumUpdateStepLength","max-step-length");
  command.AddOptionField("MaximumUpdateStepLength","floatval",MetaCommand::FLOAT,true,"2.0");

  command.SetOption("UpdateRule","a",false,"Type of update rule. 0: s <- s o exp(u) (diffeomorphic), 1: s <- s + u (additive, ITK basic), 2: s <- s o (Id+u) (compositive, Thirion's proposal?)");
  command.SetOptionLongTag("UpdateRule","update-rule");
  command.AddOptionField("UpdateRule","type",MetaCommand::INT,true,"0");
  command.SetOptionRange("UpdateRule","type","0","2");

  command.SetOption("GradienType","t",false,"Type of gradient used for computing the demons force. 0 is symmetrized, 1 is fixed image, 2 is warped moving image, 3 is mapped moving image");
  command.SetOptionLongTag("GradienType","gradient-type");
  command.AddOptionField("GradienType","type",MetaCommand::INT,true,"0");
  command.SetOptionRange("GradienType","type","0","3");

  command.SetOption("UseHistogramMatching","e",false,"Use histogram matching prior to registration (e.g. for different MR scanners)");
  command.SetOptionLongTag("UseHistogramMatching","use-histogram-matching");
  command.AddOptionField("UseHistogramMatching","boolval",MetaCommand::FLAG,false);

  command.SetOption("AlgorithmVerbosity","d",false,"Algorithm verbosity (debug level)");
  command.SetOptionLongTag("AlgorithmVerbosity","verbose");
  command.AddOptionField("AlgorithmVerbosity","intval",MetaCommand::INT,false,"1");
  command.SetOptionRange("AlgorithmVerbosity","intval","0","100"); */



  // Actually parse the command line
  if (!command.Parse(argc,argv))
    {
    exit( EXIT_FAILURE );
    }



  // Store the parsed information into a struct

  args.movingImageFile = command.GetValueAsString("MovingImageFile","filename");
  args.inputFieldFile = command.GetValueAsString("InputFieldFile","filename");
  args.outputImageFile = command.GetValueAsString("OutputImageFile","filename");

  args.useNearestInterpolation = command.GetValueAsBool("UseNearestInterpolation","boolval");

  /* args.sigmaDef = command.GetValueAsFloat("DeformationFieldSigma","floatval");
  args.sigmaUp = command.GetValueAsFloat("UpdateFieldSigma","floatval");
  args.maxStepLength = command.GetValueAsFloat("MaximumUpdateStepLength","floatval");
  args.updateRule = command.GetValueAsInt("UpdateRule","type");
  args.gradientType = command.GetValueAsInt("GradientType","type");
  args.useHistogramMatching = command.GetValueAsBool("UseHistogramMatching","boolval"); */


}


//  The following section of code implements a Command observer
//  that will monitor the evolution of the registration process.
//
template <class TPixel=float, unsigned int VImageDimension=3>
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate                         Self;
  typedef  itk::Command                                   Superclass;
  typedef  itk::SmartPointer<Self>                        Pointer;

  typedef itk::Image< TPixel, VImageDimension >           InternalImageType;
  typedef itk::Vector< TPixel, VImageDimension >          VectorPixelType;
  typedef itk::Image<  VectorPixelType, VImageDimension > DeformationFieldType;

  typedef itk::DiffeomorphicDemonsRegistrationFilter<
    InternalImageType,
    InternalImageType,
    DeformationFieldType>                                DiffeomorphicDemonsRegistrationFilterType;

  typedef itk::FastSymmetricForcesDemonsRegistrationFilter<
     InternalImageType,
     InternalImageType,
     DeformationFieldType>                                FastSymmetricForcesDemonsRegistrationFilterType;

  typedef itk::MultiResolutionPDEDeformableRegistration<
     InternalImageType, InternalImageType,
     DeformationFieldType, TPixel >                       MultiResRegistrationFilterType;

  typedef itk::DisplacementFieldJacobianDeterminantFilter<
     DeformationFieldType, TPixel>                        JacobianFilterType;

  typedef itk::MinimumMaximumImageCalculator<
     InternalImageType>                                   MinMaxFilterType;

  typedef itk::WarpHarmonicEnergyCalculator<
     DeformationFieldType>                                HarmonicEnergyCalculatorType;

  typedef itk::VectorCentralDifferenceImageFunction<
     DeformationFieldType>                                WarpGradientCalculatorType;

  typedef typename WarpGradientCalculatorType::OutputType WarpGradientType;

  itkNewMacro( Self );

  void SetTrueField(const DeformationFieldType * truefield)
    {
    m_TrueField = truefield;

    m_TrueWarpGradientCalculator = WarpGradientCalculatorType::New();
    m_TrueWarpGradientCalculator->SetInputImage( m_TrueField );

    m_CompWarpGradientCalculator =  WarpGradientCalculatorType::New();
    }

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
    if( !(itk::IterationEvent().CheckEvent( &event )) )
      {
      return;
      }

    typename DeformationFieldType::ConstPointer deffield = 0;
    unsigned int iter = -1;
    double metricbefore = -1.0;

    if ( const DiffeomorphicDemonsRegistrationFilterType * dfilter =
         dynamic_cast< const DiffeomorphicDemonsRegistrationFilterType * >( object ) )
      {
      iter = dfilter->GetElapsedIterations() - 1;
      metricbefore = dfilter->GetMetric();
      deffield = const_cast<DiffeomorphicDemonsRegistrationFilterType *>(
        dfilter)->GetDeformationField();
      }
    else if ( const FastSymmetricForcesDemonsRegistrationFilterType * ffilter =
              dynamic_cast< const FastSymmetricForcesDemonsRegistrationFilterType * >( object ) )
      {
      iter = ffilter->GetElapsedIterations() - 1;
      metricbefore = ffilter->GetMetric();
      deffield = const_cast<FastSymmetricForcesDemonsRegistrationFilterType *>(
        ffilter)->GetDeformationField();
      }
    else if ( const MultiResRegistrationFilterType * multiresfilter =
              dynamic_cast< const MultiResRegistrationFilterType * >( object ) )
      {
      std::cout<<"Finished Multi-resolution iteration :"<<multiresfilter->GetCurrentLevel()-1<<std::endl;
      std::cout<<"=============================="<<std::endl<<std::endl;
      }
    else
      {
      return;
      }

    if (deffield)
      {
      std::cout<<iter<<": MSE "<<metricbefore<<" - ";

      double fieldDist = -1.0;
      double fieldGradDist = -1.0;
      double tmp;
      if (m_TrueField)
        {
        typedef itk::ImageRegionConstIteratorWithIndex<DeformationFieldType>
           FieldIteratorType;
        FieldIteratorType currIter(
           deffield, deffield->GetLargestPossibleRegion() );
        FieldIteratorType trueIter(
           m_TrueField, deffield->GetLargestPossibleRegion() );

        m_CompWarpGradientCalculator->SetInputImage( deffield );

        fieldDist = 0.0;
        fieldGradDist = 0.0;
        for ( currIter.GoToBegin(), trueIter.GoToBegin();
              ! currIter.IsAtEnd(); ++currIter, ++trueIter )
          {
          fieldDist += (currIter.Value() - trueIter.Value()).GetSquaredNorm();

          // No need to add Id matrix here as we do a substraction
          tmp = (
             ( m_CompWarpGradientCalculator->EvaluateAtIndex(currIter.GetIndex())
               -m_TrueWarpGradientCalculator->EvaluateAtIndex(trueIter.GetIndex())
                ).GetVnlMatrix() ).frobenius_norm();
          fieldGradDist += tmp*tmp;
          }
        fieldDist = sqrt( fieldDist/ (double)(
                             deffield->GetLargestPossibleRegion().GetNumberOfPixels()) );
        fieldGradDist = sqrt( fieldGradDist/ (double)(
                                 deffield->GetLargestPossibleRegion().GetNumberOfPixels()) );

        std::cout<<"d(.,true) "<<fieldDist<<" - ";
        std::cout<<"d(.,Jac(true)) "<<fieldGradDist<<" - ";
        }

      m_HarmonicEnergyCalculator->SetImage( deffield );
      m_HarmonicEnergyCalculator->Compute();
      const double harmonicEnergy
         = m_HarmonicEnergyCalculator->GetHarmonicEnergy();
      std::cout<<"harmo. "<<harmonicEnergy<<" - ";


      m_JacobianFilter->SetInput( deffield );
      m_JacobianFilter->UpdateLargestPossibleRegion();


      const unsigned int numPix = m_JacobianFilter->
         GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels();

      TPixel* pix_start = m_JacobianFilter->GetOutput()->GetBufferPointer();
      TPixel* pix_end = pix_start + numPix;

      TPixel* jac_ptr;

      // Get percentage of det(Jac) below 0
      unsigned int jacBelowZero(0u);
      for (jac_ptr=pix_start; jac_ptr!=pix_end; ++jac_ptr)
        {
        if ( *jac_ptr<=0.0 ) ++jacBelowZero;
        }
      const double jacBelowZeroPrc = static_cast<double>(jacBelowZero)
         / static_cast<double>(numPix);


      // Get min an max jac
      const double minJac = *(std::min_element (pix_start, pix_end));
      const double maxJac = *(std::max_element (pix_start, pix_end));

      // Get some quantiles
      // We don't need the jacobian image
      // we can modify/sort it in place
      jac_ptr = pix_start + static_cast<unsigned int>(0.002*numPix);
      std::nth_element(pix_start, jac_ptr, pix_end);
      const double Q002 = *jac_ptr;

      jac_ptr = pix_start + static_cast<unsigned int>(0.01*numPix);
      std::nth_element(pix_start, jac_ptr, pix_end);
      const double Q01 = *jac_ptr;

      jac_ptr = pix_start + static_cast<unsigned int>(0.99*numPix);
      std::nth_element(pix_start, jac_ptr, pix_end);
      const double Q99 = *jac_ptr;

      jac_ptr = pix_start + static_cast<unsigned int>(0.998*numPix);
      std::nth_element(pix_start, jac_ptr, pix_end);
      const double Q998 = *jac_ptr;


      std::cout<<"max|Jac| "<<maxJac<<" - "
               <<"min|Jac| "<<minJac<<" - "
               <<"ratio(|Jac|<=0) "<<jacBelowZeroPrc<<std::endl;

      if (this->m_Fid.is_open())
        {
        if (! m_headerwritten)
          {
          this->m_Fid<<"Iteration"
                     <<", MSE before"
                     <<", Harmonic energy"
                     <<", min|Jac|"
                     <<", 0.2% |Jac|"
                     <<", 01% |Jac|"
                     <<", 99% |Jac|"
                     <<", 99.8% |Jac|"
                     <<", max|Jac|"
                     <<", ratio(|Jac|<=0)";

          if (m_TrueField)
            {
            this->m_Fid<<", dist(warp,true warp)"
                       <<", dist(Jac,true Jac)";
            }

          this->m_Fid<<std::endl;

          m_headerwritten = true;
          }

        this->m_Fid<<iter
                   <<", "<<metricbefore
                   <<", "<<harmonicEnergy
                   <<", "<<minJac
                   <<", "<<Q002
                   <<", "<<Q01
                   <<", "<<Q99
                   <<", "<<Q998
                   <<", "<<maxJac
                   <<", "<<jacBelowZeroPrc;

        if (m_TrueField)
          {
          this->m_Fid<<", "<<fieldDist
                     <<", "<<fieldGradDist;
          }

        this->m_Fid<<std::endl;
        }
      }
    }

protected:
  CommandIterationUpdate() :
     m_Fid( "metricvalues.csv" ),
     m_headerwritten(false)
    {
    m_JacobianFilter = JacobianFilterType::New();
    m_JacobianFilter->SetUseImageSpacing( true );
    m_JacobianFilter->ReleaseDataFlagOn();

    m_Minmaxfilter = MinMaxFilterType::New();

    m_HarmonicEnergyCalculator = HarmonicEnergyCalculatorType::New();

    m_TrueField = 0;
    m_TrueWarpGradientCalculator = 0;
    m_CompWarpGradientCalculator = 0;
    };

  ~CommandIterationUpdate()
    {
    this->m_Fid.close();
    }

private:
  std::ofstream m_Fid;
  bool m_headerwritten;
  typename JacobianFilterType::Pointer m_JacobianFilter;
  typename MinMaxFilterType::Pointer m_Minmaxfilter;
  typename HarmonicEnergyCalculatorType::Pointer m_HarmonicEnergyCalculator;
  typename DeformationFieldType::ConstPointer m_TrueField;
  typename WarpGradientCalculatorType::Pointer m_TrueWarpGradientCalculator;
  typename WarpGradientCalculatorType::Pointer m_CompWarpGradientCalculator;
};


template <unsigned int Dimension>
void DemonsRegistrationFunction( arguments args )
{
  // Declare the types of the images (float or double only)
  typedef float                                  PixelType;
  typedef itk::Image< PixelType, Dimension >     ImageType;

  typedef itk::Vector< PixelType, Dimension >    VectorPixelType;
  typedef typename itk::Image
     < VectorPixelType, Dimension >              DeformationFieldType;


  // Images we use
  typename ImageType::Pointer fixedImage = 0;
  typename ImageType::Pointer movingImage = 0;
  typename DeformationFieldType::Pointer inputDefField = 0;


  // Set up the file readers
  typedef itk::ImageFileReader< ImageType >            FixedImageReaderType;
  typedef itk::ImageFileReader< ImageType >            MovingImageReaderType;
  typedef itk::ImageFileReader< DeformationFieldType > FieldReaderType;
  typedef itk::TransformFileReader                     TransformReaderType;

  {//for mem allocations

  /* typename FixedImageReaderType::Pointer fixedImageReader
     = FixedImageReaderType::New(); */
  typename MovingImageReaderType::Pointer movingImageReader
     = MovingImageReaderType::New();

  /* fixedImageReader->SetFileName( args.fixedImageFile.c_str() ); */
  movingImageReader->SetFileName( args.movingImageFile.c_str() );


  // Update the reader
  try
    {
    /* fixedImageReader->Update(); */
    movingImageReader->Update();
    }
  catch( itk::ExceptionObject& err )
    {
    std::cout << "Could not read one of the input images." << std::endl;
    std::cout << err << std::endl;
    exit( EXIT_FAILURE );
    }

  if ( ! args.inputFieldFile.empty() )
    {
    // Set up the file readers
    typename FieldReaderType::Pointer fieldReader = FieldReaderType::New();
    fieldReader->SetFileName(  args.inputFieldFile.c_str() );

    // Update the reader
    try
      {
      fieldReader->Update();
      }
    catch( itk::ExceptionObject& err )
      {
      std::cout << "Could not read the input field." << std::endl;
      std::cout << err << std::endl;
      exit( EXIT_FAILURE );
      }

    inputDefField = fieldReader->GetOutput();
    inputDefField->DisconnectPipeline();
    }




	/* fixedImage = fixedImageReader->GetOutput();
	fixedImage->DisconnectPipeline(); */
	movingImage = movingImageReader->GetOutput();
	movingImage->DisconnectPipeline();


  }//end for mem allocations


  // Here we connect deformation field to the input directly
  typename DeformationFieldType::Pointer defField = inputDefField;

  // warp the result
  typedef itk::WarpImageFilter
     < ImageType, ImageType, DeformationFieldType >  WarperType;
  typename WarperType::Pointer warper = WarperType::New();
  warper->SetInput( movingImage );
  warper->SetOutputSpacing( movingImage->GetSpacing() );
  warper->SetOutputOrigin( movingImage->GetOrigin() );
  warper->SetOutputDirection( movingImage->GetDirection() );
  warper->SetDeformationField( defField );

	if ( args.useNearestInterpolation )
	{
		typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double >  NearestInterpolatorType;
		typename NearestInterpolatorType::Pointer interpolator = NearestInterpolatorType::New();
		warper->SetInterpolator(interpolator);
	}


  // Write warped image out to file
  typedef float                                OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::CastImageFilter
     < ImageType, OutputImageType >                CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  typename WriterType::Pointer      writer =  WriterType::New();
  typename CastFilterType::Pointer  caster =  CastFilterType::New();
  writer->SetFileName( args.outputImageFile.c_str() );
  caster->SetInput( warper->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->SetUseCompression( true );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject& err )
    {
    std::cout << "Unexpected error." << std::endl;
    std::cout << err << std::endl;
    exit( EXIT_FAILURE );
    }



}


int main( int argc, char *argv[] )
{
  struct arguments args;
  parseOpts (argc, argv, args);

  std::cout<<"Warping:"<<std::endl;
  std::cout<<args<<std::endl<<std::endl;

  // FIXME uncomment for debug only
  // itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);

  // Get the image dimension
  itk::ImageIOBase::Pointer imageIO;
  try
    {
    imageIO = itk::ImageIOFactory::CreateImageIO(
       args.movingImageFile.c_str(), itk::ImageIOFactory::ReadMode);
    if ( imageIO )
      {
      imageIO->SetFileName(args.movingImageFile.c_str());
      imageIO->ReadImageInformation();
      }
    else
      {
      std::cout << "Could not read the fixed image information." << std::endl;
      exit( EXIT_FAILURE );
      }
    }
  catch( itk::ExceptionObject& err )
    {
    std::cout << "Could not read the fixed image information." << std::endl;
    std::cout << err << std::endl;
    exit( EXIT_FAILURE );
    }

  switch ( imageIO->GetNumberOfDimensions() )
  {
  case 2:
    DemonsRegistrationFunction<2>(args);
    break;
  case 3:
    DemonsRegistrationFunction<3>(args);
    break;
  default:
    std::cout << "Unsuported dimension" << std::endl;
    exit( EXIT_FAILURE );
  }

 return EXIT_SUCCESS;
}
