// Simcenter STAR-CCM+ macro: star_macro_run.java
// Written by Simcenter STAR-CCM+ 17.04.007
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;
import star.vis.*;
import star.cadmodeler.*;
import star.flow.*;
import star.meshing.*;

public class star_macro_run extends StarMacro {

  public void execute() {
    execute0();
  }

  private void execute0() {

    Simulation simulation_0 = 
      getActiveSimulation();

    Scene scene_0 = 
      simulation_0.getSceneManager().createScene("3D-CAD View");

    scene_0.initializeAndWait();

    CadModel cadModel_0 = 
      ((CadModel) simulation_0.get(SolidModelManager.class).getObject("3D-CAD Model 1"));

    simulation_0.get(SolidModelManager.class).editCadModel(cadModel_0, scene_0);

    scene_0.openInteractive();

    scene_0.setAdvancedRenderingEnabled(false);

    SceneUpdate sceneUpdate_0 = 
      scene_0.getSceneUpdate();

    HardcopyProperties hardcopyProperties_0 = 
      sceneUpdate_0.getHardcopyProperties();

    hardcopyProperties_0.setCurrentResolutionWidth(25);

    hardcopyProperties_0.setCurrentResolutionHeight(25);

    hardcopyProperties_0.setCurrentResolutionWidth(1778);

    hardcopyProperties_0.setCurrentResolutionHeight(860);

    scene_0.resetCamera();

    LabCoordinateSystem labCoordinateSystem_0 = 
      simulation_0.getCoordinateSystemManager().getLabCoordinateSystem();

    cadModel_0.getFeatureManager().create3DSketches_2("C:\\Users\\kizan\\OneDrive\\Documents\\__MSU\\Fall 2022\\EGEN 436\\Final Project\\Airfoils\\CSV_naca\\NACA2117.csv", labCoordinateSystem_0, false, true);

    cadModel_0.update();

    scene_0.resetCamera();

    Sketch3D sketch3D_0 = 
      ((Sketch3D) cadModel_0.getFeature("Sketch3D 1"));

    ExtrusionMerge extrusionMerge_0 = 
      cadModel_0.getFeatureManager().createExtrusionMerge(sketch3D_0);

    extrusionMerge_0.setAutoPreview(true);

    cadModel_0.allowMakingPartDirty(false);

    scene_0.setTransparencyOverrideMode(SceneTransparencyOverride.MAKE_SCENE_TRANSPARENT);

    extrusionMerge_0.setDirectionOption(2);

    extrusionMerge_0.setExtrudedBodyTypeOption(0);

    Units units_0 = 
      ((Units) simulation_0.getUnitsManager().getObject("m"));

    extrusionMerge_0.getDistance().setValueAndUnits(0.1, units_0);

    extrusionMerge_0.getDistanceAsymmetric().setValueAndUnits(0.1, units_0);

    extrusionMerge_0.getOffsetDistance().setValueAndUnits(0.1, units_0);

    extrusionMerge_0.setDistanceOption(0);

    extrusionMerge_0.setCoordinateSystemOption(1);

    Units units_1 = 
      ((Units) simulation_0.getUnitsManager().getObject("deg"));

    extrusionMerge_0.getDraftAngle().setValueAndUnits(10.0, units_1);

    extrusionMerge_0.setDraftOption(0);

    extrusionMerge_0.setImportedCoordinateSystem(labCoordinateSystem_0);

    extrusionMerge_0.getDirectionAxis().setCoordinateSystem(labCoordinateSystem_0);

    extrusionMerge_0.getDirectionAxis().setUnits0(units_0);

    extrusionMerge_0.getDirectionAxis().setUnits1(units_0);

    extrusionMerge_0.getDirectionAxis().setUnits2(units_0);

    extrusionMerge_0.getDirectionAxis().setDefinition("");

    extrusionMerge_0.getDirectionAxis().setValue(new DoubleVector(new double[] {0.0, 0.0, 1.0}));

    extrusionMerge_0.setFace(null);

    extrusionMerge_0.setBody(null);

    extrusionMerge_0.setFeatureInputType(0);

    extrusionMerge_0.setInputFeatureEdges(new NeoObjectVector(new Object[] {}));

    extrusionMerge_0.setSketch(sketch3D_0);

    extrusionMerge_0.setInteractingBodies(new NeoObjectVector(new Object[] {}));

    extrusionMerge_0.setInteractingBodiesBodyGroups(new NeoObjectVector(new Object[] {}));

    extrusionMerge_0.setInteractingBodiesCadFilters(new NeoObjectVector(new Object[] {}));

    extrusionMerge_0.setInteractingSelectedBodies(false);

    extrusionMerge_0.setPostOption(0);

    extrusionMerge_0.setExtrusionOption(0);

    extrusionMerge_0.setIsBodyGroupCreation(false);

    cadModel_0.getFeatureManager().markDependentNotUptodate(extrusionMerge_0);

    extrusionMerge_0.markFeatureForEdit();

    cadModel_0.allowMakingPartDirty(true);

    cadModel_0.getFeatureManager().execute(extrusionMerge_0);

    scene_0.setTransparencyOverrideMode(SceneTransparencyOverride.USE_DISPLAYER_PROPERTY);

    SplineSketchPrimitive3D splineSketchPrimitive3D_0 = 
      ((SplineSketchPrimitive3D) sketch3D_0.getSketchPrimitive3D("Spline 1"));

    star.cadmodeler.Body cadmodelerBody_0 = 
      ((star.cadmodeler.Body) extrusionMerge_0.getBody(splineSketchPrimitive3D_0));

    cadmodelerBody_0.getUnNamedFacesDefaultAttributeName();

    cadmodelerBody_0.getUnNamedFacesDefaultAttributeName();

    cadmodelerBody_0.setPresentationName("wing_body");

    cadmodelerBody_0.setUnNamedFacesDefaultAttributeName("wing_face");

    simulation_0.get(SolidModelManager.class).endEditCadModel(cadModel_0);

    simulation_0.getSceneManager().deleteScenes(new NeoObjectVector(new Object[] {scene_0}));

    SolidModelPart solidModelPart_0 = 
      ((SolidModelPart) simulation_0.get(SimulationPartManager.class).getPart("domain"));

    simulation_0.get(SimulationPartManager.class).updateParts(new NeoObjectVector(new Object[] {solidModelPart_0}));

    cadModel_0.createParts(new NeoObjectVector(new Object[] {cadmodelerBody_0}), new NeoObjectVector(new Object[] {}), true, false, 1, false, false, 3, "SharpEdges", 30.0, 3, true, 1.0E-5, false);

    SubtractPartsOperation subtractPartsOperation_0 = 
      ((SubtractPartsOperation) simulation_0.get(MeshOperationManager.class).getObject("Subtract"));

    subtractPartsOperation_0.getInputGeometryObjects().setQuery(null);

    SolidModelPart solidModelPart_1 = 
      ((SolidModelPart) simulation_0.get(SimulationPartManager.class).getPart("wing_body"));

    subtractPartsOperation_0.getInputGeometryObjects().setObjects(solidModelPart_0, solidModelPart_1);

    Region region_0 = 
      simulation_0.getRegionManager().getRegion("Region");

    Boundary boundary_0 = 
      region_0.getBoundaryManager().getBoundary("Subtract.domain.inlet");

    FlowDirectionProfile flowDirectionProfile_0 = 
      boundary_0.getValues().get(FlowDirectionProfile.class);

    Units units_2 = 
      ((Units) simulation_0.getUnitsManager().getObject(""));

    flowDirectionProfile_0.getMethod(ConstantVectorProfileMethod.class).getQuantity().setComponentsAndUnits(0.9946683524878391, 0.10312549907335175, 0.0, units_2);

    VelocityMagnitudeProfile velocityMagnitudeProfile_0 = 
      boundary_0.getValues().get(VelocityMagnitudeProfile.class);

    Units units_3 = 
      ((Units) simulation_0.getUnitsManager().getObject("m/s"));

    velocityMagnitudeProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits(0.03229616232751263, units_3);

    PhysicsContinuum physicsContinuum_0 = 
      ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Physics 1"));

    VelocityProfile velocityProfile_0 = 
      physicsContinuum_0.getInitialConditions().get(VelocityProfile.class);

    velocityProfile_0.getMethod(ConstantVectorProfileMethod.class).getQuantity().setComponentsAndUnits(0.032123970573986804, 0.0033305578581787217, 0.0, units_3);

    subtractPartsOperation_0.execute();

    PrepareFor2dOperation prepareFor2dOperation_0 = 
      ((PrepareFor2dOperation) simulation_0.get(MeshOperationManager.class).getObject("Badge for 2D Meshing"));

    prepareFor2dOperation_0.execute();

    AutoMeshOperation2d autoMeshOperation2d_0 = 
      ((AutoMeshOperation2d) simulation_0.get(MeshOperationManager.class).getObject("Automated Mesh (2D)"));

    autoMeshOperation2d_0.execute();

    XyzInternalTable xyzInternalTable_0 = 
      ((XyzInternalTable) simulation_0.getTableManager().getTable("XYZ Internal Table"));

    TableUpdate tableUpdate_0 = 
      xyzInternalTable_0.getTableUpdate();

    tableUpdate_0.setBaseFilename("naca0503, all");

    XyzInternalTable xyzInternalTable_1 = 
      ((XyzInternalTable) simulation_0.getTableManager().getTable("XYZ Internal Table Rect"));

    TableUpdate tableUpdate_1 = 
      xyzInternalTable_1.getTableUpdate();

    tableUpdate_1.setBaseFilename("naca0503, rect");

    Solution solution_0 = 
      simulation_0.getSolution();

    solution_0.clearSolution(Solution.Clear.History, Solution.Clear.Fields, Solution.Clear.LagrangianDem);

    solution_0.initializeSolution();

    ResidualPlot residualPlot_0 = 
      ((ResidualPlot) simulation_0.getPlotManager().getPlot("Residuals"));

    PlotUpdate plotUpdate_0 = 
      residualPlot_0.getPlotUpdate();

    plotUpdate_0.setAnimationFilenameBase("naca2117, AoA=5.919179260104077, Re=2061.555330235035, res");

    residualPlot_0.openInteractive();

    simulation_0.getSimulationIterator().run();

    HardcopyProperties hardcopyProperties_1 = 
      plotUpdate_0.getHardcopyProperties();

    hardcopyProperties_1.setCurrentResolutionWidth(25);

    hardcopyProperties_1.setCurrentResolutionHeight(25);

    hardcopyProperties_1.setCurrentResolutionWidth(1778);

    hardcopyProperties_1.setCurrentResolutionHeight(860);

    simulation_0.get(SimulationPartManager.class).removeParts(new NeoObjectVector(new Object[] {solidModelPart_1}));

    Scene scene_1 = 
      simulation_0.getSceneManager().createScene("3D-CAD View");

    scene_1.initializeAndWait();

    simulation_0.get(SolidModelManager.class).editCadModel(cadModel_0, scene_1);

    scene_1.openInteractive();

    scene_1.setAdvancedRenderingEnabled(false);

    SceneUpdate sceneUpdate_1 = 
      scene_1.getSceneUpdate();

    HardcopyProperties hardcopyProperties_2 = 
      sceneUpdate_1.getHardcopyProperties();

    hardcopyProperties_2.setCurrentResolutionWidth(25);

    hardcopyProperties_2.setCurrentResolutionHeight(25);

    hardcopyProperties_1.setCurrentResolutionWidth(1780);

    hardcopyProperties_1.setCurrentResolutionHeight(861);

    hardcopyProperties_2.setCurrentResolutionWidth(1778);

    hardcopyProperties_2.setCurrentResolutionHeight(860);

    scene_1.resetCamera();

    cadModel_0.getFeatureManager().delete(new NeoObjectVector(new Object[] {extrusionMerge_0}), false, false, true);

    cadModel_0.getFeatureManager().delete(new NeoObjectVector(new Object[] {sketch3D_0}), false, false, true);

    simulation_0.get(SolidModelManager.class).endEditCadModel(cadModel_0);

    simulation_0.getSceneManager().deleteScenes(new NeoObjectVector(new Object[] {scene_1}));

    hardcopyProperties_1.setCurrentResolutionWidth(1778);

    hardcopyProperties_1.setCurrentResolutionHeight(860);
  }
}