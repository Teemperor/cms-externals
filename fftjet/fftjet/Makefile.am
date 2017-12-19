lib_LTLIBRARIES = libfftjet.la

libfftjet_la_SOURCES = AbsKernel2d.cc Kernels.cc quartic_lib.cc \
         Kernel2dFactory.cc rectangleOverlap.cc StatAccumulator.cc Peak.cc \
         AbsSymmetricKernel.cc ProfileKernel.cc LogProfileKernel.cc \
         PhiKernels.cc InterpolatedKernel.cc PythiaKernel_30_100_v0.cc \
         PeakFinder.cc LinearInterpolator1d.cc LinearInterpolator2d.cc \
         cdfCellNumber.cc CompositeKernel.cc MultiKernelConvolver.cc \
         EquidistantSequence.cc InterpolatedKernel3d.cc SpecialFunctions.cc \
         DXGlyphInfo.cc PeakEtaPhiDistance.cc DiscreteGauss1d.cc \
         DiscreteGauss2d.cc AbsKernel1d.cc Kernels1d.cc Kernel1dFactory.cc \
         PeakEtaDependentDistance.cc

libfftjet_la_CPPFLAGS = -I@top_srcdir@/ $(DEPS_CFLAGS)

libfftjet_la_LIBADD = $(DEPS_LIBS)

includedir = ${prefix}/include/fftjet

include_HEADERS = AbsClusteringTree.hh \
    AbsClusteringTree.icc \
    AbsConvolverBase.hh \
    AbsDistanceCalculator.hh \
    AbsFFTEngine.hh \
    AbsFrequencyKernel.hh \
    AbsFrequencyKernel1d.hh \
    AbsJetCorrector.hh \
    AbsKernel1d.hh \
    AbsKernel1d.icc \
    AbsKernel2d.hh \
    AbsKernel2d.icc \
    AbsKernelConvolver.hh \
    AbsKernelConvolver.icc \
    AbsMembershipFunction.hh \
    AbsOpenDXTreeFormatter.hh \
    AbsOpenDXTreeFormatter.icc \
    AbsPatternRecognitionAlg.hh \
    AbsPeakSelector.hh \
    AbsPhiKernel.hh \
    AbsRecombinationAlg.hh \
    AbsRecombinationAlg.icc \
    AbsSequentialConvolver.hh \
    AbsSequentialConvolver.icc \
    AbsVectorRecombinationAlg.hh \
    AbsScalableKernel.hh \
    AbsScalableKernel1d.hh \
    AbsScalablePhiKernel.hh \
    AbsSymmetricKernel.hh \
    AbsTreeFormatter.hh \
    binaryIO.hh \
    cdfCellNumber.hh \
    ClusteringSequencer.hh \
    ClusteringSequencer.icc \
    ClusteringTreeSparsifier.hh \
    ClusteringTreeSparsifier.icc \
    CompositeKernel.hh \
    ConstScaleReconstruction.hh \
    ConstScaleReconstruction.icc \
    convertPeakCoords.hh \
    CUFFTFloatEngine.hh \
    DiscreteGauss1d.hh \
    DiscreteGauss2d.hh \
    DXGlyphInfo.hh \
    EquidistantSequence.hh \
    EtaPhiEtInterpolator.hh \
    EtaPhiEtInterpolator.icc \
    EtCentroidRecombinationAlg.hh \
    EtCentroidVectorRecombinationAlg.hh \
    EtSumRecombinationAlg.hh \
    EtSumVectorRecombinationAlg.hh \
    FasterEtCentroidRecombinationAlg.hh \
    FasterEtSumRecombinationAlg.hh \
    FasterKernelRecombinationAlg.hh \
    FasterKernelRecombinationAlg.icc \
    FFTWDoubleEngine.hh \
    FFTWDoubleEngine.icc \
    FFTWEngine.hh \
    FFTWEngine.icc \
    FFTWFloatEngine.hh \
    FFTWFloatEngine.icc \
    FrequencyKernelConvolver.hh \
    FrequencyKernelConvolver.icc \
    FrequencySequentialConvolver.hh \
    FrequencySequentialConvolver.icc \
    GaussianNoiseMembershipFcn.hh \
    Grid2d.hh \
    Grid2d.icc \
    InterpolatedKernel3d.hh \
    InterpolatedKernel3d.icc \
    InterpolatedKernel.hh \
    InterpolatedMembershipFcn.hh \
    InterpolatedMembershipFcn.icc \
    interpolate.hh \
    invertJetResponse.hh \
    invertJetResponse.icc \
    JetMagnitudeMapper2d.hh \
    JetMagnitudeMapper2d.icc \
    JetMagnitudeMapper.hh \
    JetMagnitudeMapper.icc \
    JetProperty.hh \
    Kernel1dFactory.hh \
    Kernel2dFactory.hh \
    KernelConvolver.hh \
    KernelConvolver.icc \
    KernelData.hh \
    KernelData.icc \
    KernelRecombinationAlg.hh \
    KernelRecombinationAlg.icc \
    KernelVectorRecombinationAlg.hh \
    KernelVectorRecombinationAlg.icc \
    Kernels.hh \
    Kernels1d.hh \
    LinearInterpolator1d.hh \
    LinearInterpolator1d.icc \
    LinearInterpolator2d.hh \
    LinearInterpolator2d.icc \
    LogProfileKernel.hh \
    MagneticSmearingKernel.hh \
    MagneticSmearingKernel.icc \
    Matrix.hh \
    Matrix.icc \
    MultiKernelConvolver.hh \
    MultiKernelConvolver.icc \
    OpenDXPeakTree.hh \
    Peak.hh \
    PeakEtaPhiDistance.hh \
    PeakEtaDependentDistance.hh \
    PeakFinder.hh \
    PeakFinder.icc \
    peakEtLifetime.hh \
    peakEtLifetime.icc \
    PeakSelectors.hh \
    PeakSelectors.icc \
    PhiKernels.hh \
    ProfileKernel.hh \
    ProximityClusteringTree.hh \
    ProximityClusteringTree.icc \
    PythiaKernel_30_100_v0.hh \
    quartic_lib.hh \
    RealFrequencyKernel.hh \
    RealFrequencyKernel1d.hh \
    RecombinationAlgFactory.hh \
    RecombinationAlgFactory.icc \
    RecombinedJet.hh \
    RecombinedJet.icc \
    rectangleOverlap.hh \
    ScaleSpaceKernel.hh \
    scanFrequencyKernel.hh \
    SequentialConvolver.hh \
    SequentialConvolver.icc \
    SimpleFunctors.hh \
    SimplePredicates.hh \
    SmallVector.hh \
    SmallVector.icc \
    SparseClusteringTree.hh \
    SparseClusteringTree.icc \
    SpecialFunctions.hh \
    StatAccumulator.hh \
    VectorRecombinationAlgFactory.hh \
    VectorRecombinationAlgFactory.icc

EXTRA_DIST = pythia_kernel_30_100_v0.scc rectangleOverlap.scc \
             CUFFTFloatEngine.cc
