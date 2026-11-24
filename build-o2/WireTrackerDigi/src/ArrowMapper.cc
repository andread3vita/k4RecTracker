// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "podio/utilities/ArrowFrameConverter.h"
#include "podio/utilities/ArrowTypeRegistry.h"
#include "podio/utilities/ArrowConverterRegistry.h"
#include "podio/CollectionBufferFactory.h"
#include <arrow/api.h>

#include "podio/utilities/ArrowUtils.h"

#include <map>
#include <memory>
#include <stdexcept>
#include <string>

#include "podio/LinkCollection.h"

#include "extension/DriftChamberDigiLocalCollection.h"
#include "extension/DriftChamberDigiCollection.h"
#include "extension/MCRecoDriftChamberDigiAssociationCollection.h"
#include "extension/DriftChamberDigiV2Collection.h"
#include "extension/MCRecoDriftChamberDigiV2AssociationCollection.h"
#include "extension/TrackCollection.h"
#include "extension/SenseWireHitCollection.h"
#include "extension/SenseWireHitSimTrackerHitLinkCollection.h"

namespace extension::arrow_io {

namespace {
using podio::arrow_utils::checkStatus;



static const std::map<std::string, std::shared_ptr<arrow::DataType>>& typeMap() {
  static const std::map<std::string, std::shared_ptr<arrow::DataType>> instance = {
    {
      "extension::DriftChamberDigiLocal",
      arrow::list(arrow::struct_({
          arrow::field("cellID", arrow::uint64(), false),
          arrow::field("distanceToWire", arrow::float32(), false),
          arrow::field("zPositionAlongWire", arrow::float32(), false),
          arrow::field("time", arrow::float32(), false),
          arrow::field("eDep", arrow::float32(), false),
          arrow::field("eDepError", arrow::float32(), false),
          arrow::field("clusterCount", arrow::uint32(), false),
      }))
    },
    {
      "extension::DriftChamberDigi",
      arrow::list(arrow::struct_({
          arrow::field("cellID", arrow::uint64(), false),
          arrow::field("leftPosition", arrow::struct_({arrow::field("x", arrow::float64(), false), arrow::field("y", arrow::float64(), false), arrow::field("z", arrow::float64(), false)}), false),
          arrow::field("rightPosition", arrow::struct_({arrow::field("x", arrow::float64(), false), arrow::field("y", arrow::float64(), false), arrow::field("z", arrow::float64(), false)}), false),
          arrow::field("time", arrow::float32(), false),
          arrow::field("eDep", arrow::float32(), false),
          arrow::field("eDepError", arrow::float32(), false),
          arrow::field("clusterCount", arrow::uint32(), false),
          arrow::field("type", arrow::int32(), false),
          arrow::field("quality", arrow::int32(), false),
          arrow::field("position", arrow::struct_({arrow::field("x", arrow::float64(), false), arrow::field("y", arrow::float64(), false), arrow::field("z", arrow::float64(), false)}), false),
      }))
    },
    {
      "extension::MCRecoDriftChamberDigiAssociation",
      arrow::list(arrow::struct_({
          arrow::field("weight", arrow::float32(), false),
          arrow::field("digi", podio::objectRefType(), false),
          arrow::field("sim", podio::objectRefType(), false),
      }))
    },
    {
      "extension::DriftChamberDigiV2",
      arrow::list(arrow::struct_({
          arrow::field("cellID", arrow::uint64(), false),
          arrow::field("type", arrow::int32(), false),
          arrow::field("quality", arrow::int32(), false),
          arrow::field("time", arrow::float32(), false),
          arrow::field("eDep", arrow::float32(), false),
          arrow::field("eDepError", arrow::float32(), false),
          arrow::field("position", arrow::struct_({arrow::field("x", arrow::float64(), false), arrow::field("y", arrow::float64(), false), arrow::field("z", arrow::float64(), false)}), false),
          arrow::field("directionSW", arrow::struct_({arrow::field("x", arrow::float64(), false), arrow::field("y", arrow::float64(), false), arrow::field("z", arrow::float64(), false)}), false),
          arrow::field("distanceToWire", arrow::float32(), false),
          arrow::field("nCluster", arrow::uint32(), false),
          arrow::field("nElectrons", arrow::list(arrow::uint16()), false),
      }))
    },
    {
      "extension::MCRecoDriftChamberDigiV2Association",
      arrow::list(arrow::struct_({
          arrow::field("weight", arrow::float32(), false),
          arrow::field("digi", podio::objectRefType(), false),
          arrow::field("sim", podio::objectRefType(), false),
      }))
    },
    {
      "extension::Track",
      arrow::list(arrow::struct_({
          arrow::field("type", arrow::int32(), false),
          arrow::field("chi2", arrow::float32(), false),
          arrow::field("ndf", arrow::int32(), false),
          arrow::field("dEdx", arrow::float32(), false),
          arrow::field("dEdxError", arrow::float32(), false),
          arrow::field("radiusOfInnermostHit", arrow::float32(), false),
          arrow::field("subdetectorHitNumbers", arrow::list(arrow::int32()), false),
          arrow::field("trackStates", arrow::list(arrow::struct_({arrow::field("location", arrow::int32(), false), arrow::field("D0", arrow::float32(), false), arrow::field("phi", arrow::float32(), false), arrow::field("omega", arrow::float32(), false), arrow::field("Z0", arrow::float32(), false), arrow::field("tanLambda", arrow::float32(), false), arrow::field("time", arrow::float32(), false), arrow::field("referencePoint", arrow::struct_({arrow::field("x", arrow::float32(), false), arrow::field("y", arrow::float32(), false), arrow::field("z", arrow::float32(), false)}), false), arrow::field("covMatrix", arrow::struct_({arrow::field("values", arrow::fixed_size_list(arrow::float32(), 21), false)}), false)})), false),
          arrow::field("dxQuantities", arrow::list(arrow::struct_({arrow::field("type", arrow::int32(), false), arrow::field("value", arrow::float32(), false), arrow::field("error", arrow::float32(), false)})), false),
          arrow::field("trackerHits", arrow::list(podio::objectRefType()), false),
          arrow::field("tracks", arrow::list(podio::objectRefType()), false),
      }))
    },
    {
      "extension::SenseWireHit",
      arrow::list(arrow::struct_({
          arrow::field("cellID", arrow::uint64(), false),
          arrow::field("type", arrow::int32(), false),
          arrow::field("quality", arrow::int32(), false),
          arrow::field("time", arrow::float32(), false),
          arrow::field("eDep", arrow::float32(), false),
          arrow::field("eDepError", arrow::float32(), false),
          arrow::field("wireStereoAngle", arrow::float32(), false),
          arrow::field("wireAzimuthalAngle", arrow::float32(), false),
          arrow::field("position", arrow::struct_({arrow::field("x", arrow::float64(), false), arrow::field("y", arrow::float64(), false), arrow::field("z", arrow::float64(), false)}), false),
          arrow::field("positionAlongWireError", arrow::float64(), false),
          arrow::field("distanceToWire", arrow::float32(), false),
          arrow::field("distanceToWireError", arrow::float32(), false),
          arrow::field("nElectrons", arrow::list(arrow::uint16()), false),
      }))
    },
    {
      "extension::SenseWireHitSimTrackerHitLink",
      arrow::list(arrow::struct_({
          arrow::field("weight", arrow::float32(), false),
          arrow::field("from", podio::objectRefType(), false),
          arrow::field("to", podio::objectRefType(), false),
      }))
    },
  };
  return instance;
}

bool registerTypes() {
  auto& reg = podio::ArrowTypeRegistry::mutInstance();
  for (const auto& [typeName, type] : typeMap()) {
    reg.registerType(typeName, type);
  }
  return true;
}

static const bool typesRegistered = registerTypes();

bool registerConverters() {
  auto& convRegistry = podio::ArrowConverterRegistry::mutInstance();
  convRegistry.registerConverter("extension::DriftChamberDigiLocal", [](const podio::CollectionBase* coll) {
    const auto* concreteColl = static_cast<const extension::DriftChamberDigiLocalCollection*>(coll);
    auto type = podio::ArrowTypeRegistry::instance().getType("extension::DriftChamberDigiLocal");
    std::unique_ptr<arrow::ArrayBuilder> builder;
    auto status = arrow::MakeBuilder(arrow::default_memory_pool(), type, &builder);
    if (!status.ok()) {
      throw std::runtime_error("Failed to create builder for type 'extension::DriftChamberDigiLocal': " + status.ToString());
    }
    auto* collectionBuilder = static_cast<arrow::ListBuilder*>(builder.get());
    auto* objectBuilder = static_cast<arrow::StructBuilder*>(collectionBuilder->value_builder());
    auto* builder_cellID = static_cast<arrow::UInt64Builder*>(objectBuilder->child(0));
    auto* builder_distanceToWire = static_cast<arrow::FloatBuilder*>(objectBuilder->child(1));
    auto* builder_zPositionAlongWire = static_cast<arrow::FloatBuilder*>(objectBuilder->child(2));
    auto* builder_time = static_cast<arrow::FloatBuilder*>(objectBuilder->child(3));
    auto* builder_eDep = static_cast<arrow::FloatBuilder*>(objectBuilder->child(4));
    auto* builder_eDepError = static_cast<arrow::FloatBuilder*>(objectBuilder->child(5));
    auto* builder_clusterCount = static_cast<arrow::UInt32Builder*>(objectBuilder->child(6));
    checkStatus(collectionBuilder->Append(), "Failed to append to collectionBuilder");
    for (const auto obj: *concreteColl) {
      checkStatus(builder_cellID->Append(obj.getCellID()), "Failed to append member cellID");
      checkStatus(builder_distanceToWire->Append(obj.getDistanceToWire()), "Failed to append member distanceToWire");
      checkStatus(builder_zPositionAlongWire->Append(obj.getZPositionAlongWire()), "Failed to append member zPositionAlongWire");
      checkStatus(builder_time->Append(obj.getTime()), "Failed to append member time");
      checkStatus(builder_eDep->Append(obj.getEDep()), "Failed to append member eDep");
      checkStatus(builder_eDepError->Append(obj.getEDepError()), "Failed to append member eDepError");
      checkStatus(builder_clusterCount->Append(obj.getClusterCount()), "Failed to append member clusterCount");
      checkStatus(objectBuilder->Append(), "Failed to append to objectBuilder");
    }
    std::shared_ptr<arrow::Array> array;
    checkStatus(collectionBuilder->Finish(&array), "Failed to finish collectionBuilder");
    return array;
  });
  convRegistry.registerConverter("extension::DriftChamberDigi", [](const podio::CollectionBase* coll) {
    const auto* concreteColl = static_cast<const extension::DriftChamberDigiCollection*>(coll);
    auto type = podio::ArrowTypeRegistry::instance().getType("extension::DriftChamberDigi");
    std::unique_ptr<arrow::ArrayBuilder> builder;
    auto status = arrow::MakeBuilder(arrow::default_memory_pool(), type, &builder);
    if (!status.ok()) {
      throw std::runtime_error("Failed to create builder for type 'extension::DriftChamberDigi': " + status.ToString());
    }
    auto* collectionBuilder = static_cast<arrow::ListBuilder*>(builder.get());
    auto* objectBuilder = static_cast<arrow::StructBuilder*>(collectionBuilder->value_builder());
    auto* builder_cellID = static_cast<arrow::UInt64Builder*>(objectBuilder->child(0));
    auto* builder_leftPosition = static_cast<arrow::StructBuilder*>(objectBuilder->child(1));
    auto* builder_leftPosition_x = static_cast<arrow::DoubleBuilder*>(builder_leftPosition->child(0));

    auto* builder_leftPosition_y = static_cast<arrow::DoubleBuilder*>(builder_leftPosition->child(1));

    auto* builder_leftPosition_z = static_cast<arrow::DoubleBuilder*>(builder_leftPosition->child(2));

    auto* builder_rightPosition = static_cast<arrow::StructBuilder*>(objectBuilder->child(2));
    auto* builder_rightPosition_x = static_cast<arrow::DoubleBuilder*>(builder_rightPosition->child(0));

    auto* builder_rightPosition_y = static_cast<arrow::DoubleBuilder*>(builder_rightPosition->child(1));

    auto* builder_rightPosition_z = static_cast<arrow::DoubleBuilder*>(builder_rightPosition->child(2));

    auto* builder_time = static_cast<arrow::FloatBuilder*>(objectBuilder->child(3));
    auto* builder_eDep = static_cast<arrow::FloatBuilder*>(objectBuilder->child(4));
    auto* builder_eDepError = static_cast<arrow::FloatBuilder*>(objectBuilder->child(5));
    auto* builder_clusterCount = static_cast<arrow::UInt32Builder*>(objectBuilder->child(6));
    auto* builder_type = static_cast<arrow::Int32Builder*>(objectBuilder->child(7));
    auto* builder_quality = static_cast<arrow::Int32Builder*>(objectBuilder->child(8));
    auto* builder_position = static_cast<arrow::StructBuilder*>(objectBuilder->child(9));
    auto* builder_position_x = static_cast<arrow::DoubleBuilder*>(builder_position->child(0));

    auto* builder_position_y = static_cast<arrow::DoubleBuilder*>(builder_position->child(1));

    auto* builder_position_z = static_cast<arrow::DoubleBuilder*>(builder_position->child(2));

    checkStatus(collectionBuilder->Append(), "Failed to append to collectionBuilder");
    for (const auto obj: *concreteColl) {
      checkStatus(builder_cellID->Append(obj.getCellID()), "Failed to append member cellID");
      checkStatus(builder_leftPosition->Append(), "Failed to append struct");
      checkStatus(builder_leftPosition_x->Append(obj.getLeftPosition().x), "Failed to append member x");

      checkStatus(builder_leftPosition_y->Append(obj.getLeftPosition().y), "Failed to append member y");

      checkStatus(builder_leftPosition_z->Append(obj.getLeftPosition().z), "Failed to append member z");

      checkStatus(builder_rightPosition->Append(), "Failed to append struct");
      checkStatus(builder_rightPosition_x->Append(obj.getRightPosition().x), "Failed to append member x");

      checkStatus(builder_rightPosition_y->Append(obj.getRightPosition().y), "Failed to append member y");

      checkStatus(builder_rightPosition_z->Append(obj.getRightPosition().z), "Failed to append member z");

      checkStatus(builder_time->Append(obj.getTime()), "Failed to append member time");
      checkStatus(builder_eDep->Append(obj.getEDep()), "Failed to append member eDep");
      checkStatus(builder_eDepError->Append(obj.getEDepError()), "Failed to append member eDepError");
      checkStatus(builder_clusterCount->Append(obj.getClusterCount()), "Failed to append member clusterCount");
      checkStatus(builder_type->Append(obj.getType()), "Failed to append member type");
      checkStatus(builder_quality->Append(obj.getQuality()), "Failed to append member quality");
      checkStatus(builder_position->Append(), "Failed to append struct");
      checkStatus(builder_position_x->Append(obj.getPosition().x), "Failed to append member x");

      checkStatus(builder_position_y->Append(obj.getPosition().y), "Failed to append member y");

      checkStatus(builder_position_z->Append(obj.getPosition().z), "Failed to append member z");

      checkStatus(objectBuilder->Append(), "Failed to append to objectBuilder");
    }
    std::shared_ptr<arrow::Array> array;
    checkStatus(collectionBuilder->Finish(&array), "Failed to finish collectionBuilder");
    return array;
  });
  convRegistry.registerConverter("extension::MCRecoDriftChamberDigiAssociation", [](const podio::CollectionBase* coll) {
    const auto* concreteColl = static_cast<const extension::MCRecoDriftChamberDigiAssociationCollection*>(coll);
    auto type = podio::ArrowTypeRegistry::instance().getType("extension::MCRecoDriftChamberDigiAssociation");
    std::unique_ptr<arrow::ArrayBuilder> builder;
    auto status = arrow::MakeBuilder(arrow::default_memory_pool(), type, &builder);
    if (!status.ok()) {
      throw std::runtime_error("Failed to create builder for type 'extension::MCRecoDriftChamberDigiAssociation': " + status.ToString());
    }
    auto* collectionBuilder = static_cast<arrow::ListBuilder*>(builder.get());
    auto* objectBuilder = static_cast<arrow::StructBuilder*>(collectionBuilder->value_builder());
    auto* builder_weight = static_cast<arrow::FloatBuilder*>(objectBuilder->child(0));
    auto* builder_digi = static_cast<arrow::StructBuilder*>(objectBuilder->child(1));
    auto* builder_digi_collId = static_cast<arrow::UInt32Builder*>(builder_digi->child(0));
    auto* builder_digi_index = static_cast<arrow::Int32Builder*>(builder_digi->child(1));
    auto* builder_sim = static_cast<arrow::StructBuilder*>(objectBuilder->child(2));
    auto* builder_sim_collId = static_cast<arrow::UInt32Builder*>(builder_sim->child(0));
    auto* builder_sim_index = static_cast<arrow::Int32Builder*>(builder_sim->child(1));
    checkStatus(collectionBuilder->Append(), "Failed to append to collectionBuilder");
    for (const auto obj: *concreteColl) {
      checkStatus(builder_weight->Append(obj.getWeight()), "Failed to append member weight");
      auto relObj_digi = obj.getDigi();
      checkStatus(builder_digi->Append(), "Failed to append relation struct");
      checkStatus(builder_digi_collId->Append(relObj_digi.getObjectID().collectionID), "Failed to append relation collectionID");
      checkStatus(builder_digi_index->Append(relObj_digi.getObjectID().index), "Failed to append relation index");
      auto relObj_sim = obj.getSim();
      checkStatus(builder_sim->Append(), "Failed to append relation struct");
      checkStatus(builder_sim_collId->Append(relObj_sim.getObjectID().collectionID), "Failed to append relation collectionID");
      checkStatus(builder_sim_index->Append(relObj_sim.getObjectID().index), "Failed to append relation index");
      checkStatus(objectBuilder->Append(), "Failed to append to objectBuilder");
    }
    std::shared_ptr<arrow::Array> array;
    checkStatus(collectionBuilder->Finish(&array), "Failed to finish collectionBuilder");
    return array;
  });
  convRegistry.registerConverter("extension::DriftChamberDigiV2", [](const podio::CollectionBase* coll) {
    const auto* concreteColl = static_cast<const extension::DriftChamberDigiV2Collection*>(coll);
    auto type = podio::ArrowTypeRegistry::instance().getType("extension::DriftChamberDigiV2");
    std::unique_ptr<arrow::ArrayBuilder> builder;
    auto status = arrow::MakeBuilder(arrow::default_memory_pool(), type, &builder);
    if (!status.ok()) {
      throw std::runtime_error("Failed to create builder for type 'extension::DriftChamberDigiV2': " + status.ToString());
    }
    auto* collectionBuilder = static_cast<arrow::ListBuilder*>(builder.get());
    auto* objectBuilder = static_cast<arrow::StructBuilder*>(collectionBuilder->value_builder());
    auto* builder_cellID = static_cast<arrow::UInt64Builder*>(objectBuilder->child(0));
    auto* builder_type = static_cast<arrow::Int32Builder*>(objectBuilder->child(1));
    auto* builder_quality = static_cast<arrow::Int32Builder*>(objectBuilder->child(2));
    auto* builder_time = static_cast<arrow::FloatBuilder*>(objectBuilder->child(3));
    auto* builder_eDep = static_cast<arrow::FloatBuilder*>(objectBuilder->child(4));
    auto* builder_eDepError = static_cast<arrow::FloatBuilder*>(objectBuilder->child(5));
    auto* builder_position = static_cast<arrow::StructBuilder*>(objectBuilder->child(6));
    auto* builder_position_x = static_cast<arrow::DoubleBuilder*>(builder_position->child(0));

    auto* builder_position_y = static_cast<arrow::DoubleBuilder*>(builder_position->child(1));

    auto* builder_position_z = static_cast<arrow::DoubleBuilder*>(builder_position->child(2));

    auto* builder_directionSW = static_cast<arrow::StructBuilder*>(objectBuilder->child(7));
    auto* builder_directionSW_x = static_cast<arrow::DoubleBuilder*>(builder_directionSW->child(0));

    auto* builder_directionSW_y = static_cast<arrow::DoubleBuilder*>(builder_directionSW->child(1));

    auto* builder_directionSW_z = static_cast<arrow::DoubleBuilder*>(builder_directionSW->child(2));

    auto* builder_distanceToWire = static_cast<arrow::FloatBuilder*>(objectBuilder->child(8));
    auto* builder_nCluster = static_cast<arrow::UInt32Builder*>(objectBuilder->child(9));
    auto* builder_nElectrons = static_cast<arrow::ListBuilder*>(objectBuilder->child(10));
    auto* builder_nElectrons_val = static_cast<arrow::UInt16Builder*>(builder_nElectrons->value_builder());
    checkStatus(collectionBuilder->Append(), "Failed to append to collectionBuilder");
    for (const auto obj: *concreteColl) {
      checkStatus(builder_cellID->Append(obj.getCellID()), "Failed to append member cellID");
      checkStatus(builder_type->Append(obj.getType()), "Failed to append member type");
      checkStatus(builder_quality->Append(obj.getQuality()), "Failed to append member quality");
      checkStatus(builder_time->Append(obj.getTime()), "Failed to append member time");
      checkStatus(builder_eDep->Append(obj.getEDep()), "Failed to append member eDep");
      checkStatus(builder_eDepError->Append(obj.getEDepError()), "Failed to append member eDepError");
      checkStatus(builder_position->Append(), "Failed to append struct");
      checkStatus(builder_position_x->Append(obj.getPosition().x), "Failed to append member x");

      checkStatus(builder_position_y->Append(obj.getPosition().y), "Failed to append member y");

      checkStatus(builder_position_z->Append(obj.getPosition().z), "Failed to append member z");

      checkStatus(builder_directionSW->Append(), "Failed to append struct");
      checkStatus(builder_directionSW_x->Append(obj.getDirectionSW().x), "Failed to append member x");

      checkStatus(builder_directionSW_y->Append(obj.getDirectionSW().y), "Failed to append member y");

      checkStatus(builder_directionSW_z->Append(obj.getDirectionSW().z), "Failed to append member z");

      checkStatus(builder_distanceToWire->Append(obj.getDistanceToWire()), "Failed to append member distanceToWire");
      checkStatus(builder_nCluster->Append(obj.getNCluster()), "Failed to append member nCluster");
      checkStatus(builder_nElectrons->Append(), "Failed to append list");
      for (size_t vec_i = 0; vec_i < obj.nElectrons_size(); ++vec_i) {
      checkStatus(builder_nElectrons_val->Append(obj.getNElectrons(vec_i)), "Failed to append vector element");
      }      checkStatus(objectBuilder->Append(), "Failed to append to objectBuilder");
    }
    std::shared_ptr<arrow::Array> array;
    checkStatus(collectionBuilder->Finish(&array), "Failed to finish collectionBuilder");
    return array;
  });
  convRegistry.registerConverter("extension::MCRecoDriftChamberDigiV2Association", [](const podio::CollectionBase* coll) {
    const auto* concreteColl = static_cast<const extension::MCRecoDriftChamberDigiV2AssociationCollection*>(coll);
    auto type = podio::ArrowTypeRegistry::instance().getType("extension::MCRecoDriftChamberDigiV2Association");
    std::unique_ptr<arrow::ArrayBuilder> builder;
    auto status = arrow::MakeBuilder(arrow::default_memory_pool(), type, &builder);
    if (!status.ok()) {
      throw std::runtime_error("Failed to create builder for type 'extension::MCRecoDriftChamberDigiV2Association': " + status.ToString());
    }
    auto* collectionBuilder = static_cast<arrow::ListBuilder*>(builder.get());
    auto* objectBuilder = static_cast<arrow::StructBuilder*>(collectionBuilder->value_builder());
    auto* builder_weight = static_cast<arrow::FloatBuilder*>(objectBuilder->child(0));
    auto* builder_digi = static_cast<arrow::StructBuilder*>(objectBuilder->child(1));
    auto* builder_digi_collId = static_cast<arrow::UInt32Builder*>(builder_digi->child(0));
    auto* builder_digi_index = static_cast<arrow::Int32Builder*>(builder_digi->child(1));
    auto* builder_sim = static_cast<arrow::StructBuilder*>(objectBuilder->child(2));
    auto* builder_sim_collId = static_cast<arrow::UInt32Builder*>(builder_sim->child(0));
    auto* builder_sim_index = static_cast<arrow::Int32Builder*>(builder_sim->child(1));
    checkStatus(collectionBuilder->Append(), "Failed to append to collectionBuilder");
    for (const auto obj: *concreteColl) {
      checkStatus(builder_weight->Append(obj.getWeight()), "Failed to append member weight");
      auto relObj_digi = obj.getDigi();
      checkStatus(builder_digi->Append(), "Failed to append relation struct");
      checkStatus(builder_digi_collId->Append(relObj_digi.getObjectID().collectionID), "Failed to append relation collectionID");
      checkStatus(builder_digi_index->Append(relObj_digi.getObjectID().index), "Failed to append relation index");
      auto relObj_sim = obj.getSim();
      checkStatus(builder_sim->Append(), "Failed to append relation struct");
      checkStatus(builder_sim_collId->Append(relObj_sim.getObjectID().collectionID), "Failed to append relation collectionID");
      checkStatus(builder_sim_index->Append(relObj_sim.getObjectID().index), "Failed to append relation index");
      checkStatus(objectBuilder->Append(), "Failed to append to objectBuilder");
    }
    std::shared_ptr<arrow::Array> array;
    checkStatus(collectionBuilder->Finish(&array), "Failed to finish collectionBuilder");
    return array;
  });
  convRegistry.registerConverter("extension::Track", [](const podio::CollectionBase* coll) {
    const auto* concreteColl = static_cast<const extension::TrackCollection*>(coll);
    auto type = podio::ArrowTypeRegistry::instance().getType("extension::Track");
    std::unique_ptr<arrow::ArrayBuilder> builder;
    auto status = arrow::MakeBuilder(arrow::default_memory_pool(), type, &builder);
    if (!status.ok()) {
      throw std::runtime_error("Failed to create builder for type 'extension::Track': " + status.ToString());
    }
    auto* collectionBuilder = static_cast<arrow::ListBuilder*>(builder.get());
    auto* objectBuilder = static_cast<arrow::StructBuilder*>(collectionBuilder->value_builder());
    auto* builder_type = static_cast<arrow::Int32Builder*>(objectBuilder->child(0));
    auto* builder_chi2 = static_cast<arrow::FloatBuilder*>(objectBuilder->child(1));
    auto* builder_ndf = static_cast<arrow::Int32Builder*>(objectBuilder->child(2));
    auto* builder_dEdx = static_cast<arrow::FloatBuilder*>(objectBuilder->child(3));
    auto* builder_dEdxError = static_cast<arrow::FloatBuilder*>(objectBuilder->child(4));
    auto* builder_radiusOfInnermostHit = static_cast<arrow::FloatBuilder*>(objectBuilder->child(5));
    auto* builder_subdetectorHitNumbers = static_cast<arrow::ListBuilder*>(objectBuilder->child(6));
    auto* builder_subdetectorHitNumbers_val = static_cast<arrow::Int32Builder*>(builder_subdetectorHitNumbers->value_builder());
    auto* builder_trackStates = static_cast<arrow::ListBuilder*>(objectBuilder->child(7));
    auto* builder_trackStates_val = static_cast<arrow::StructBuilder*>(builder_trackStates->value_builder());
    auto* builder_trackStates_val_location = static_cast<arrow::Int32Builder*>(builder_trackStates_val->child(0));

    auto* builder_trackStates_val_D0 = static_cast<arrow::FloatBuilder*>(builder_trackStates_val->child(1));

    auto* builder_trackStates_val_phi = static_cast<arrow::FloatBuilder*>(builder_trackStates_val->child(2));

    auto* builder_trackStates_val_omega = static_cast<arrow::FloatBuilder*>(builder_trackStates_val->child(3));

    auto* builder_trackStates_val_Z0 = static_cast<arrow::FloatBuilder*>(builder_trackStates_val->child(4));

    auto* builder_trackStates_val_tanLambda = static_cast<arrow::FloatBuilder*>(builder_trackStates_val->child(5));

    auto* builder_trackStates_val_time = static_cast<arrow::FloatBuilder*>(builder_trackStates_val->child(6));

    auto* builder_trackStates_val_referencePoint = static_cast<arrow::StructBuilder*>(builder_trackStates_val->child(7));
    auto* builder_trackStates_val_referencePoint_x = static_cast<arrow::FloatBuilder*>(builder_trackStates_val_referencePoint->child(0));

    auto* builder_trackStates_val_referencePoint_y = static_cast<arrow::FloatBuilder*>(builder_trackStates_val_referencePoint->child(1));

    auto* builder_trackStates_val_referencePoint_z = static_cast<arrow::FloatBuilder*>(builder_trackStates_val_referencePoint->child(2));


    auto* builder_trackStates_val_covMatrix = static_cast<arrow::StructBuilder*>(builder_trackStates_val->child(8));
    auto* builder_trackStates_val_covMatrix_values = static_cast<arrow::FixedSizeListBuilder*>(builder_trackStates_val_covMatrix->child(0));
    auto* builder_trackStates_val_covMatrix_values_val = static_cast<arrow::FloatBuilder*>(builder_trackStates_val_covMatrix_values->value_builder());


    auto* builder_dxQuantities = static_cast<arrow::ListBuilder*>(objectBuilder->child(8));
    auto* builder_dxQuantities_val = static_cast<arrow::StructBuilder*>(builder_dxQuantities->value_builder());
    auto* builder_dxQuantities_val_type = static_cast<arrow::Int32Builder*>(builder_dxQuantities_val->child(0));

    auto* builder_dxQuantities_val_value = static_cast<arrow::FloatBuilder*>(builder_dxQuantities_val->child(1));

    auto* builder_dxQuantities_val_error = static_cast<arrow::FloatBuilder*>(builder_dxQuantities_val->child(2));

    auto* builder_trackerHits = static_cast<arrow::ListBuilder*>(objectBuilder->child(9));
    auto* builder_trackerHits_struct = static_cast<arrow::StructBuilder*>(builder_trackerHits->value_builder());
    auto* builder_trackerHits_collId = static_cast<arrow::UInt32Builder*>(builder_trackerHits_struct->child(0));
    auto* builder_trackerHits_index = static_cast<arrow::Int32Builder*>(builder_trackerHits_struct->child(1));
    auto* builder_tracks = static_cast<arrow::ListBuilder*>(objectBuilder->child(10));
    auto* builder_tracks_struct = static_cast<arrow::StructBuilder*>(builder_tracks->value_builder());
    auto* builder_tracks_collId = static_cast<arrow::UInt32Builder*>(builder_tracks_struct->child(0));
    auto* builder_tracks_index = static_cast<arrow::Int32Builder*>(builder_tracks_struct->child(1));
    checkStatus(collectionBuilder->Append(), "Failed to append to collectionBuilder");
    for (const auto obj: *concreteColl) {
      checkStatus(builder_type->Append(obj.getType()), "Failed to append member type");
      checkStatus(builder_chi2->Append(obj.getChi2()), "Failed to append member chi2");
      checkStatus(builder_ndf->Append(obj.getNdf()), "Failed to append member ndf");
      checkStatus(builder_dEdx->Append(obj.getDEdx()), "Failed to append member dEdx");
      checkStatus(builder_dEdxError->Append(obj.getDEdxError()), "Failed to append member dEdxError");
      checkStatus(builder_radiusOfInnermostHit->Append(obj.getRadiusOfInnermostHit()), "Failed to append member radiusOfInnermostHit");
      checkStatus(builder_subdetectorHitNumbers->Append(), "Failed to append list");
      for (size_t vec_i = 0; vec_i < obj.subdetectorHitNumbers_size(); ++vec_i) {
      checkStatus(builder_subdetectorHitNumbers_val->Append(obj.getSubdetectorHitNumbers(vec_i)), "Failed to append vector element");
      }      checkStatus(builder_trackStates->Append(), "Failed to append list");
      for (size_t vec_i = 0; vec_i < obj.trackStates_size(); ++vec_i) {
      checkStatus(builder_trackStates_val->Append(), "Failed to append vector struct element");
        checkStatus(builder_trackStates_val_location->Append(obj.getTrackStates(vec_i).location), "Failed to append member location");

        checkStatus(builder_trackStates_val_D0->Append(obj.getTrackStates(vec_i).D0), "Failed to append member D0");

        checkStatus(builder_trackStates_val_phi->Append(obj.getTrackStates(vec_i).phi), "Failed to append member phi");

        checkStatus(builder_trackStates_val_omega->Append(obj.getTrackStates(vec_i).omega), "Failed to append member omega");

        checkStatus(builder_trackStates_val_Z0->Append(obj.getTrackStates(vec_i).Z0), "Failed to append member Z0");

        checkStatus(builder_trackStates_val_tanLambda->Append(obj.getTrackStates(vec_i).tanLambda), "Failed to append member tanLambda");

        checkStatus(builder_trackStates_val_time->Append(obj.getTrackStates(vec_i).time), "Failed to append member time");

        checkStatus(builder_trackStates_val_referencePoint->Append(), "Failed to append struct");
        checkStatus(builder_trackStates_val_referencePoint_x->Append(obj.getTrackStates(vec_i).referencePoint.x), "Failed to append member x");

        checkStatus(builder_trackStates_val_referencePoint_y->Append(obj.getTrackStates(vec_i).referencePoint.y), "Failed to append member y");

        checkStatus(builder_trackStates_val_referencePoint_z->Append(obj.getTrackStates(vec_i).referencePoint.z), "Failed to append member z");


        checkStatus(builder_trackStates_val_covMatrix->Append(), "Failed to append struct");
        checkStatus(builder_trackStates_val_covMatrix_values->Append(), "Failed to append fixed_size_list");
        for (size_t arr_i = 0; arr_i < 21; ++arr_i) {
        checkStatus(builder_trackStates_val_covMatrix_values_val->Append(obj.getTrackStates(vec_i).covMatrix.values[arr_i]), "Failed to append array element");
        }


      }      checkStatus(builder_dxQuantities->Append(), "Failed to append list");
      for (size_t vec_i = 0; vec_i < obj.dxQuantities_size(); ++vec_i) {
      checkStatus(builder_dxQuantities_val->Append(), "Failed to append vector struct element");
        checkStatus(builder_dxQuantities_val_type->Append(obj.getDxQuantities(vec_i).type), "Failed to append member type");

        checkStatus(builder_dxQuantities_val_value->Append(obj.getDxQuantities(vec_i).value), "Failed to append member value");

        checkStatus(builder_dxQuantities_val_error->Append(obj.getDxQuantities(vec_i).error), "Failed to append member error");

      }      checkStatus(builder_trackerHits->Append(), "Failed to append relation list");
      for (size_t rel_i = 0; rel_i < obj.trackerHits_size(); ++rel_i) {
        auto relObj = obj.getTrackerHits(rel_i);
        checkStatus(builder_trackerHits_struct->Append(), "Failed to append relation list struct");
        checkStatus(builder_trackerHits_collId->Append(relObj.getObjectID().collectionID), "Failed to append relation list collectionID");
        checkStatus(builder_trackerHits_index->Append(relObj.getObjectID().index), "Failed to append relation list index");
      }
      checkStatus(builder_tracks->Append(), "Failed to append relation list");
      for (size_t rel_i = 0; rel_i < obj.tracks_size(); ++rel_i) {
        auto relObj = obj.getTracks(rel_i);
        checkStatus(builder_tracks_struct->Append(), "Failed to append relation list struct");
        checkStatus(builder_tracks_collId->Append(relObj.getObjectID().collectionID), "Failed to append relation list collectionID");
        checkStatus(builder_tracks_index->Append(relObj.getObjectID().index), "Failed to append relation list index");
      }
      checkStatus(objectBuilder->Append(), "Failed to append to objectBuilder");
    }
    std::shared_ptr<arrow::Array> array;
    checkStatus(collectionBuilder->Finish(&array), "Failed to finish collectionBuilder");
    return array;
  });
  convRegistry.registerConverter("extension::SenseWireHit", [](const podio::CollectionBase* coll) {
    const auto* concreteColl = static_cast<const extension::SenseWireHitCollection*>(coll);
    auto type = podio::ArrowTypeRegistry::instance().getType("extension::SenseWireHit");
    std::unique_ptr<arrow::ArrayBuilder> builder;
    auto status = arrow::MakeBuilder(arrow::default_memory_pool(), type, &builder);
    if (!status.ok()) {
      throw std::runtime_error("Failed to create builder for type 'extension::SenseWireHit': " + status.ToString());
    }
    auto* collectionBuilder = static_cast<arrow::ListBuilder*>(builder.get());
    auto* objectBuilder = static_cast<arrow::StructBuilder*>(collectionBuilder->value_builder());
    auto* builder_cellID = static_cast<arrow::UInt64Builder*>(objectBuilder->child(0));
    auto* builder_type = static_cast<arrow::Int32Builder*>(objectBuilder->child(1));
    auto* builder_quality = static_cast<arrow::Int32Builder*>(objectBuilder->child(2));
    auto* builder_time = static_cast<arrow::FloatBuilder*>(objectBuilder->child(3));
    auto* builder_eDep = static_cast<arrow::FloatBuilder*>(objectBuilder->child(4));
    auto* builder_eDepError = static_cast<arrow::FloatBuilder*>(objectBuilder->child(5));
    auto* builder_wireStereoAngle = static_cast<arrow::FloatBuilder*>(objectBuilder->child(6));
    auto* builder_wireAzimuthalAngle = static_cast<arrow::FloatBuilder*>(objectBuilder->child(7));
    auto* builder_position = static_cast<arrow::StructBuilder*>(objectBuilder->child(8));
    auto* builder_position_x = static_cast<arrow::DoubleBuilder*>(builder_position->child(0));

    auto* builder_position_y = static_cast<arrow::DoubleBuilder*>(builder_position->child(1));

    auto* builder_position_z = static_cast<arrow::DoubleBuilder*>(builder_position->child(2));

    auto* builder_positionAlongWireError = static_cast<arrow::DoubleBuilder*>(objectBuilder->child(9));
    auto* builder_distanceToWire = static_cast<arrow::FloatBuilder*>(objectBuilder->child(10));
    auto* builder_distanceToWireError = static_cast<arrow::FloatBuilder*>(objectBuilder->child(11));
    auto* builder_nElectrons = static_cast<arrow::ListBuilder*>(objectBuilder->child(12));
    auto* builder_nElectrons_val = static_cast<arrow::UInt16Builder*>(builder_nElectrons->value_builder());
    checkStatus(collectionBuilder->Append(), "Failed to append to collectionBuilder");
    for (const auto obj: *concreteColl) {
      checkStatus(builder_cellID->Append(obj.getCellID()), "Failed to append member cellID");
      checkStatus(builder_type->Append(obj.getType()), "Failed to append member type");
      checkStatus(builder_quality->Append(obj.getQuality()), "Failed to append member quality");
      checkStatus(builder_time->Append(obj.getTime()), "Failed to append member time");
      checkStatus(builder_eDep->Append(obj.getEDep()), "Failed to append member eDep");
      checkStatus(builder_eDepError->Append(obj.getEDepError()), "Failed to append member eDepError");
      checkStatus(builder_wireStereoAngle->Append(obj.getWireStereoAngle()), "Failed to append member wireStereoAngle");
      checkStatus(builder_wireAzimuthalAngle->Append(obj.getWireAzimuthalAngle()), "Failed to append member wireAzimuthalAngle");
      checkStatus(builder_position->Append(), "Failed to append struct");
      checkStatus(builder_position_x->Append(obj.getPosition().x), "Failed to append member x");

      checkStatus(builder_position_y->Append(obj.getPosition().y), "Failed to append member y");

      checkStatus(builder_position_z->Append(obj.getPosition().z), "Failed to append member z");

      checkStatus(builder_positionAlongWireError->Append(obj.getPositionAlongWireError()), "Failed to append member positionAlongWireError");
      checkStatus(builder_distanceToWire->Append(obj.getDistanceToWire()), "Failed to append member distanceToWire");
      checkStatus(builder_distanceToWireError->Append(obj.getDistanceToWireError()), "Failed to append member distanceToWireError");
      checkStatus(builder_nElectrons->Append(), "Failed to append list");
      for (size_t vec_i = 0; vec_i < obj.nElectrons_size(); ++vec_i) {
      checkStatus(builder_nElectrons_val->Append(obj.getNElectrons(vec_i)), "Failed to append vector element");
      }      checkStatus(objectBuilder->Append(), "Failed to append to objectBuilder");
    }
    std::shared_ptr<arrow::Array> array;
    checkStatus(collectionBuilder->Finish(&array), "Failed to finish collectionBuilder");
    return array;
  });
  convRegistry.registerConverter("extension::SenseWireHitSimTrackerHitLink", [](const podio::CollectionBase* coll) {
    const auto* concreteColl = static_cast<const extension::SenseWireHitSimTrackerHitLinkCollection*>(coll);
    auto type = podio::ArrowTypeRegistry::instance().getType("extension::SenseWireHitSimTrackerHitLink");
    std::unique_ptr<arrow::ArrayBuilder> builder;
    auto status = arrow::MakeBuilder(arrow::default_memory_pool(), type, &builder);
    if (!status.ok()) {
      throw std::runtime_error("Failed to create builder for type 'extension::SenseWireHitSimTrackerHitLink': " + status.ToString());
    }
    auto* collectionBuilder = static_cast<arrow::ListBuilder*>(builder.get());
    auto* objectBuilder = static_cast<arrow::StructBuilder*>(collectionBuilder->value_builder());
    auto* builder_weight = static_cast<arrow::FloatBuilder*>(objectBuilder->child(0));
    auto* builder_from = static_cast<arrow::StructBuilder*>(objectBuilder->child(1));
    auto* builder_from_collId = static_cast<arrow::UInt32Builder*>(builder_from->child(0));
    auto* builder_from_index = static_cast<arrow::Int32Builder*>(builder_from->child(1));
    auto* builder_to = static_cast<arrow::StructBuilder*>(objectBuilder->child(2));
    auto* builder_to_collId = static_cast<arrow::UInt32Builder*>(builder_to->child(0));
    auto* builder_to_index = static_cast<arrow::Int32Builder*>(builder_to->child(1));
    checkStatus(collectionBuilder->Append(), "Failed to append to collectionBuilder");
    for (const auto obj: *concreteColl) {
      checkStatus(builder_weight->Append(obj.getWeight()), "Failed to append member weight");
      auto relObj_from = obj.getFrom();
      checkStatus(builder_from->Append(), "Failed to append relation struct");
      checkStatus(builder_from_collId->Append(relObj_from.getObjectID().collectionID), "Failed to append relation collectionID");
      checkStatus(builder_from_index->Append(relObj_from.getObjectID().index), "Failed to append relation index");
      auto relObj_to = obj.getTo();
      checkStatus(builder_to->Append(), "Failed to append relation struct");
      checkStatus(builder_to_collId->Append(relObj_to.getObjectID().collectionID), "Failed to append relation collectionID");
      checkStatus(builder_to_index->Append(relObj_to.getObjectID().index), "Failed to append relation index");
      checkStatus(objectBuilder->Append(), "Failed to append to objectBuilder");
    }
    std::shared_ptr<arrow::Array> array;
    checkStatus(collectionBuilder->Finish(&array), "Failed to finish collectionBuilder");
    return array;
  });
  return true;
}

static const bool convertersRegistered = registerConverters();

bool registerBufferReaders() {
  auto& convRegistry = podio::ArrowConverterRegistry::mutInstance();
  convRegistry.registerReader("extension::DriftChamberDigiLocal", [](std::shared_ptr<arrow::Array> array, int64_t rowIndex, bool isSubset, podio::SchemaVersionT version) -> std::optional<podio::CollectionReadBuffers> {
    auto buffers = podio::CollectionBufferFactory::instance().createBuffers("extension::DriftChamberDigiLocalCollection", version, isSubset);
    if (!buffers) {
      return std::nullopt;
    }

    auto collection_frames_list = std::static_pointer_cast<arrow::ListArray>(array);
    auto frame_collection_slice = collection_frames_list->value_slice(rowIndex);
    auto raw_data_struct_array = std::static_pointer_cast<arrow::StructArray>(frame_collection_slice);

    if (!buffers->data) {
      auto collId_arr = std::static_pointer_cast<arrow::UInt32Array>(raw_data_struct_array->field(0));
      auto index_arr = std::static_pointer_cast<arrow::Int32Array>(raw_data_struct_array->field(1));
      auto* refVec = buffers->references->at(0).get();
      size_t collection_size = raw_data_struct_array->length();
      refVec->reserve(collection_size);
      for (size_t i = 0; i < collection_size; ++i) {
        refVec->push_back(podio::ObjectID{index_arr->Value(i), collId_arr->Value(i)});
      }
      return buffers;
    }

    auto* dataVec = buffers->dataAsVector<extension::DriftChamberDigiLocalData>();
    size_t collection_size = raw_data_struct_array->length();
    dataVec->resize(collection_size);    auto arr_cellID = std::static_pointer_cast<arrow::UInt64Array>(raw_data_struct_array->GetFieldByName("cellID"));
    auto arr_distanceToWire = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("distanceToWire"));
    auto arr_zPositionAlongWire = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("zPositionAlongWire"));
    auto arr_time = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("time"));
    auto arr_eDep = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("eDep"));
    auto arr_eDepError = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("eDepError"));
    auto arr_clusterCount = std::static_pointer_cast<arrow::UInt32Array>(raw_data_struct_array->GetFieldByName("clusterCount"));

    for (size_t i = 0; i < collection_size; ++i) {      (*dataVec)[i].cellID = arr_cellID->Value(i);
      (*dataVec)[i].distanceToWire = arr_distanceToWire->Value(i);
      (*dataVec)[i].zPositionAlongWire = arr_zPositionAlongWire->Value(i);
      (*dataVec)[i].time = arr_time->Value(i);
      (*dataVec)[i].eDep = arr_eDep->Value(i);
      (*dataVec)[i].eDepError = arr_eDepError->Value(i);
      (*dataVec)[i].clusterCount = arr_clusterCount->Value(i);
    }
    return buffers;
  });
  convRegistry.registerReader("extension::DriftChamberDigi", [](std::shared_ptr<arrow::Array> array, int64_t rowIndex, bool isSubset, podio::SchemaVersionT version) -> std::optional<podio::CollectionReadBuffers> {
    auto buffers = podio::CollectionBufferFactory::instance().createBuffers("extension::DriftChamberDigiCollection", version, isSubset);
    if (!buffers) {
      return std::nullopt;
    }

    auto collection_frames_list = std::static_pointer_cast<arrow::ListArray>(array);
    auto frame_collection_slice = collection_frames_list->value_slice(rowIndex);
    auto raw_data_struct_array = std::static_pointer_cast<arrow::StructArray>(frame_collection_slice);

    if (!buffers->data) {
      auto collId_arr = std::static_pointer_cast<arrow::UInt32Array>(raw_data_struct_array->field(0));
      auto index_arr = std::static_pointer_cast<arrow::Int32Array>(raw_data_struct_array->field(1));
      auto* refVec = buffers->references->at(0).get();
      size_t collection_size = raw_data_struct_array->length();
      refVec->reserve(collection_size);
      for (size_t i = 0; i < collection_size; ++i) {
        refVec->push_back(podio::ObjectID{index_arr->Value(i), collId_arr->Value(i)});
      }
      return buffers;
    }

    auto* dataVec = buffers->dataAsVector<extension::DriftChamberDigiData>();
    size_t collection_size = raw_data_struct_array->length();
    dataVec->resize(collection_size);    auto arr_cellID = std::static_pointer_cast<arrow::UInt64Array>(raw_data_struct_array->GetFieldByName("cellID"));
    auto arr_leftPosition = std::static_pointer_cast<arrow::StructArray>(raw_data_struct_array->GetFieldByName("leftPosition"));
    auto arr_leftPosition_x = std::static_pointer_cast<arrow::DoubleArray>(arr_leftPosition->GetFieldByName("x"));

    auto arr_leftPosition_y = std::static_pointer_cast<arrow::DoubleArray>(arr_leftPosition->GetFieldByName("y"));

    auto arr_leftPosition_z = std::static_pointer_cast<arrow::DoubleArray>(arr_leftPosition->GetFieldByName("z"));

    auto arr_rightPosition = std::static_pointer_cast<arrow::StructArray>(raw_data_struct_array->GetFieldByName("rightPosition"));
    auto arr_rightPosition_x = std::static_pointer_cast<arrow::DoubleArray>(arr_rightPosition->GetFieldByName("x"));

    auto arr_rightPosition_y = std::static_pointer_cast<arrow::DoubleArray>(arr_rightPosition->GetFieldByName("y"));

    auto arr_rightPosition_z = std::static_pointer_cast<arrow::DoubleArray>(arr_rightPosition->GetFieldByName("z"));

    auto arr_time = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("time"));
    auto arr_eDep = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("eDep"));
    auto arr_eDepError = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("eDepError"));
    auto arr_clusterCount = std::static_pointer_cast<arrow::UInt32Array>(raw_data_struct_array->GetFieldByName("clusterCount"));
    auto arr_type = std::static_pointer_cast<arrow::Int32Array>(raw_data_struct_array->GetFieldByName("type"));
    auto arr_quality = std::static_pointer_cast<arrow::Int32Array>(raw_data_struct_array->GetFieldByName("quality"));
    auto arr_position = std::static_pointer_cast<arrow::StructArray>(raw_data_struct_array->GetFieldByName("position"));
    auto arr_position_x = std::static_pointer_cast<arrow::DoubleArray>(arr_position->GetFieldByName("x"));

    auto arr_position_y = std::static_pointer_cast<arrow::DoubleArray>(arr_position->GetFieldByName("y"));

    auto arr_position_z = std::static_pointer_cast<arrow::DoubleArray>(arr_position->GetFieldByName("z"));


    for (size_t i = 0; i < collection_size; ++i) {      (*dataVec)[i].cellID = arr_cellID->Value(i);
      (*dataVec)[i].leftPosition.x = arr_leftPosition_x->Value(i);

      (*dataVec)[i].leftPosition.y = arr_leftPosition_y->Value(i);

      (*dataVec)[i].leftPosition.z = arr_leftPosition_z->Value(i);

      (*dataVec)[i].rightPosition.x = arr_rightPosition_x->Value(i);

      (*dataVec)[i].rightPosition.y = arr_rightPosition_y->Value(i);

      (*dataVec)[i].rightPosition.z = arr_rightPosition_z->Value(i);

      (*dataVec)[i].time = arr_time->Value(i);
      (*dataVec)[i].eDep = arr_eDep->Value(i);
      (*dataVec)[i].eDepError = arr_eDepError->Value(i);
      (*dataVec)[i].clusterCount = arr_clusterCount->Value(i);
      (*dataVec)[i].type = arr_type->Value(i);
      (*dataVec)[i].quality = arr_quality->Value(i);
      (*dataVec)[i].position.x = arr_position_x->Value(i);

      (*dataVec)[i].position.y = arr_position_y->Value(i);

      (*dataVec)[i].position.z = arr_position_z->Value(i);

    }
    return buffers;
  });
  convRegistry.registerReader("extension::MCRecoDriftChamberDigiAssociation", [](std::shared_ptr<arrow::Array> array, int64_t rowIndex, bool isSubset, podio::SchemaVersionT version) -> std::optional<podio::CollectionReadBuffers> {
    auto buffers = podio::CollectionBufferFactory::instance().createBuffers("extension::MCRecoDriftChamberDigiAssociationCollection", version, isSubset);
    if (!buffers) {
      return std::nullopt;
    }

    auto collection_frames_list = std::static_pointer_cast<arrow::ListArray>(array);
    auto frame_collection_slice = collection_frames_list->value_slice(rowIndex);
    auto raw_data_struct_array = std::static_pointer_cast<arrow::StructArray>(frame_collection_slice);

    if (!buffers->data) {
      auto collId_arr = std::static_pointer_cast<arrow::UInt32Array>(raw_data_struct_array->field(0));
      auto index_arr = std::static_pointer_cast<arrow::Int32Array>(raw_data_struct_array->field(1));
      auto* refVec = buffers->references->at(0).get();
      size_t collection_size = raw_data_struct_array->length();
      refVec->reserve(collection_size);
      for (size_t i = 0; i < collection_size; ++i) {
        refVec->push_back(podio::ObjectID{index_arr->Value(i), collId_arr->Value(i)});
      }
      return buffers;
    }

    auto* dataVec = buffers->dataAsVector<extension::MCRecoDriftChamberDigiAssociationData>();
    size_t collection_size = raw_data_struct_array->length();
    dataVec->resize(collection_size);    auto arr_weight = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("weight"));
    auto rel_struct_arr_digi = std::static_pointer_cast<arrow::StructArray>(raw_data_struct_array->GetFieldByName("digi"));
    auto rel_collId_arr_digi = std::static_pointer_cast<arrow::UInt32Array>(rel_struct_arr_digi->field(0));
    auto rel_index_arr_digi = std::static_pointer_cast<arrow::Int32Array>(rel_struct_arr_digi->field(1));
    auto* refVec_digi = buffers->references->at(0).get();
    refVec_digi->reserve(collection_size);    auto rel_struct_arr_sim = std::static_pointer_cast<arrow::StructArray>(raw_data_struct_array->GetFieldByName("sim"));
    auto rel_collId_arr_sim = std::static_pointer_cast<arrow::UInt32Array>(rel_struct_arr_sim->field(0));
    auto rel_index_arr_sim = std::static_pointer_cast<arrow::Int32Array>(rel_struct_arr_sim->field(1));
    auto* refVec_sim = buffers->references->at(1).get();
    refVec_sim->reserve(collection_size);
    for (size_t i = 0; i < collection_size; ++i) {      (*dataVec)[i].weight = arr_weight->Value(i);
      {
        refVec_digi->push_back(podio::ObjectID{rel_index_arr_digi->Value(i), rel_collId_arr_digi->Value(i)});
      }      {
        refVec_sim->push_back(podio::ObjectID{rel_index_arr_sim->Value(i), rel_collId_arr_sim->Value(i)});
      }    }
    return buffers;
  });
  convRegistry.registerReader("extension::DriftChamberDigiV2", [](std::shared_ptr<arrow::Array> array, int64_t rowIndex, bool isSubset, podio::SchemaVersionT version) -> std::optional<podio::CollectionReadBuffers> {
    auto buffers = podio::CollectionBufferFactory::instance().createBuffers("extension::DriftChamberDigiV2Collection", version, isSubset);
    if (!buffers) {
      return std::nullopt;
    }

    auto collection_frames_list = std::static_pointer_cast<arrow::ListArray>(array);
    auto frame_collection_slice = collection_frames_list->value_slice(rowIndex);
    auto raw_data_struct_array = std::static_pointer_cast<arrow::StructArray>(frame_collection_slice);

    if (!buffers->data) {
      auto collId_arr = std::static_pointer_cast<arrow::UInt32Array>(raw_data_struct_array->field(0));
      auto index_arr = std::static_pointer_cast<arrow::Int32Array>(raw_data_struct_array->field(1));
      auto* refVec = buffers->references->at(0).get();
      size_t collection_size = raw_data_struct_array->length();
      refVec->reserve(collection_size);
      for (size_t i = 0; i < collection_size; ++i) {
        refVec->push_back(podio::ObjectID{index_arr->Value(i), collId_arr->Value(i)});
      }
      return buffers;
    }

    auto* dataVec = buffers->dataAsVector<extension::DriftChamberDigiV2Data>();
    size_t collection_size = raw_data_struct_array->length();
    dataVec->resize(collection_size);    auto arr_cellID = std::static_pointer_cast<arrow::UInt64Array>(raw_data_struct_array->GetFieldByName("cellID"));
    auto arr_type = std::static_pointer_cast<arrow::Int32Array>(raw_data_struct_array->GetFieldByName("type"));
    auto arr_quality = std::static_pointer_cast<arrow::Int32Array>(raw_data_struct_array->GetFieldByName("quality"));
    auto arr_time = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("time"));
    auto arr_eDep = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("eDep"));
    auto arr_eDepError = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("eDepError"));
    auto arr_position = std::static_pointer_cast<arrow::StructArray>(raw_data_struct_array->GetFieldByName("position"));
    auto arr_position_x = std::static_pointer_cast<arrow::DoubleArray>(arr_position->GetFieldByName("x"));

    auto arr_position_y = std::static_pointer_cast<arrow::DoubleArray>(arr_position->GetFieldByName("y"));

    auto arr_position_z = std::static_pointer_cast<arrow::DoubleArray>(arr_position->GetFieldByName("z"));

    auto arr_directionSW = std::static_pointer_cast<arrow::StructArray>(raw_data_struct_array->GetFieldByName("directionSW"));
    auto arr_directionSW_x = std::static_pointer_cast<arrow::DoubleArray>(arr_directionSW->GetFieldByName("x"));

    auto arr_directionSW_y = std::static_pointer_cast<arrow::DoubleArray>(arr_directionSW->GetFieldByName("y"));

    auto arr_directionSW_z = std::static_pointer_cast<arrow::DoubleArray>(arr_directionSW->GetFieldByName("z"));

    auto arr_distanceToWire = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("distanceToWire"));
    auto arr_nCluster = std::static_pointer_cast<arrow::UInt32Array>(raw_data_struct_array->GetFieldByName("nCluster"));
    auto arr_nElectrons = std::static_pointer_cast<arrow::ListArray>(raw_data_struct_array->GetFieldByName("nElectrons"));    auto arr_nElectrons_val = std::static_pointer_cast<arrow::UInt16Array>(arr_nElectrons->values());    auto* vecMemberBuffer_nElectrons = podio::CollectionReadBuffers::asVector<std::uint16_t>(buffers->vectorMembers->at(0).second);
    vecMemberBuffer_nElectrons->reserve(arr_nElectrons->values()->length());
    for (size_t i = 0; i < collection_size; ++i) {      (*dataVec)[i].cellID = arr_cellID->Value(i);
      (*dataVec)[i].type = arr_type->Value(i);
      (*dataVec)[i].quality = arr_quality->Value(i);
      (*dataVec)[i].time = arr_time->Value(i);
      (*dataVec)[i].eDep = arr_eDep->Value(i);
      (*dataVec)[i].eDepError = arr_eDepError->Value(i);
      (*dataVec)[i].position.x = arr_position_x->Value(i);

      (*dataVec)[i].position.y = arr_position_y->Value(i);

      (*dataVec)[i].position.z = arr_position_z->Value(i);

      (*dataVec)[i].directionSW.x = arr_directionSW_x->Value(i);

      (*dataVec)[i].directionSW.y = arr_directionSW_y->Value(i);

      (*dataVec)[i].directionSW.z = arr_directionSW_z->Value(i);

      (*dataVec)[i].distanceToWire = arr_distanceToWire->Value(i);
      (*dataVec)[i].nCluster = arr_nCluster->Value(i);
      {
        size_t count_begin = vecMemberBuffer_nElectrons->size();
        // Get start and end offsets of the vector/list element in the flat values array
        int32_t start_offset = arr_nElectrons->value_offset(i);
        int32_t end_offset = arr_nElectrons->value_offset(i + 1);
        for (int32_t j = start_offset; j < end_offset; ++j) {
vecMemberBuffer_nElectrons->push_back(arr_nElectrons_val->Value(j));
}
        (*dataVec)[i].nElectrons_begin = count_begin;
        (*dataVec)[i].nElectrons_end = vecMemberBuffer_nElectrons->size();
      }    }
    return buffers;
  });
  convRegistry.registerReader("extension::MCRecoDriftChamberDigiV2Association", [](std::shared_ptr<arrow::Array> array, int64_t rowIndex, bool isSubset, podio::SchemaVersionT version) -> std::optional<podio::CollectionReadBuffers> {
    auto buffers = podio::CollectionBufferFactory::instance().createBuffers("extension::MCRecoDriftChamberDigiV2AssociationCollection", version, isSubset);
    if (!buffers) {
      return std::nullopt;
    }

    auto collection_frames_list = std::static_pointer_cast<arrow::ListArray>(array);
    auto frame_collection_slice = collection_frames_list->value_slice(rowIndex);
    auto raw_data_struct_array = std::static_pointer_cast<arrow::StructArray>(frame_collection_slice);

    if (!buffers->data) {
      auto collId_arr = std::static_pointer_cast<arrow::UInt32Array>(raw_data_struct_array->field(0));
      auto index_arr = std::static_pointer_cast<arrow::Int32Array>(raw_data_struct_array->field(1));
      auto* refVec = buffers->references->at(0).get();
      size_t collection_size = raw_data_struct_array->length();
      refVec->reserve(collection_size);
      for (size_t i = 0; i < collection_size; ++i) {
        refVec->push_back(podio::ObjectID{index_arr->Value(i), collId_arr->Value(i)});
      }
      return buffers;
    }

    auto* dataVec = buffers->dataAsVector<extension::MCRecoDriftChamberDigiV2AssociationData>();
    size_t collection_size = raw_data_struct_array->length();
    dataVec->resize(collection_size);    auto arr_weight = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("weight"));
    auto rel_struct_arr_digi = std::static_pointer_cast<arrow::StructArray>(raw_data_struct_array->GetFieldByName("digi"));
    auto rel_collId_arr_digi = std::static_pointer_cast<arrow::UInt32Array>(rel_struct_arr_digi->field(0));
    auto rel_index_arr_digi = std::static_pointer_cast<arrow::Int32Array>(rel_struct_arr_digi->field(1));
    auto* refVec_digi = buffers->references->at(0).get();
    refVec_digi->reserve(collection_size);    auto rel_struct_arr_sim = std::static_pointer_cast<arrow::StructArray>(raw_data_struct_array->GetFieldByName("sim"));
    auto rel_collId_arr_sim = std::static_pointer_cast<arrow::UInt32Array>(rel_struct_arr_sim->field(0));
    auto rel_index_arr_sim = std::static_pointer_cast<arrow::Int32Array>(rel_struct_arr_sim->field(1));
    auto* refVec_sim = buffers->references->at(1).get();
    refVec_sim->reserve(collection_size);
    for (size_t i = 0; i < collection_size; ++i) {      (*dataVec)[i].weight = arr_weight->Value(i);
      {
        refVec_digi->push_back(podio::ObjectID{rel_index_arr_digi->Value(i), rel_collId_arr_digi->Value(i)});
      }      {
        refVec_sim->push_back(podio::ObjectID{rel_index_arr_sim->Value(i), rel_collId_arr_sim->Value(i)});
      }    }
    return buffers;
  });
  convRegistry.registerReader("extension::Track", [](std::shared_ptr<arrow::Array> array, int64_t rowIndex, bool isSubset, podio::SchemaVersionT version) -> std::optional<podio::CollectionReadBuffers> {
    auto buffers = podio::CollectionBufferFactory::instance().createBuffers("extension::TrackCollection", version, isSubset);
    if (!buffers) {
      return std::nullopt;
    }

    auto collection_frames_list = std::static_pointer_cast<arrow::ListArray>(array);
    auto frame_collection_slice = collection_frames_list->value_slice(rowIndex);
    auto raw_data_struct_array = std::static_pointer_cast<arrow::StructArray>(frame_collection_slice);

    if (!buffers->data) {
      auto collId_arr = std::static_pointer_cast<arrow::UInt32Array>(raw_data_struct_array->field(0));
      auto index_arr = std::static_pointer_cast<arrow::Int32Array>(raw_data_struct_array->field(1));
      auto* refVec = buffers->references->at(0).get();
      size_t collection_size = raw_data_struct_array->length();
      refVec->reserve(collection_size);
      for (size_t i = 0; i < collection_size; ++i) {
        refVec->push_back(podio::ObjectID{index_arr->Value(i), collId_arr->Value(i)});
      }
      return buffers;
    }

    auto* dataVec = buffers->dataAsVector<extension::TrackData>();
    size_t collection_size = raw_data_struct_array->length();
    dataVec->resize(collection_size);    auto arr_type = std::static_pointer_cast<arrow::Int32Array>(raw_data_struct_array->GetFieldByName("type"));
    auto arr_chi2 = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("chi2"));
    auto arr_ndf = std::static_pointer_cast<arrow::Int32Array>(raw_data_struct_array->GetFieldByName("ndf"));
    auto arr_dEdx = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("dEdx"));
    auto arr_dEdxError = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("dEdxError"));
    auto arr_radiusOfInnermostHit = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("radiusOfInnermostHit"));
    auto arr_subdetectorHitNumbers = std::static_pointer_cast<arrow::ListArray>(raw_data_struct_array->GetFieldByName("subdetectorHitNumbers"));    auto arr_subdetectorHitNumbers_val = std::static_pointer_cast<arrow::Int32Array>(arr_subdetectorHitNumbers->values());    auto* vecMemberBuffer_subdetectorHitNumbers = podio::CollectionReadBuffers::asVector<std::int32_t>(buffers->vectorMembers->at(0).second);
    vecMemberBuffer_subdetectorHitNumbers->reserve(arr_subdetectorHitNumbers->values()->length());    auto arr_trackStates = std::static_pointer_cast<arrow::ListArray>(raw_data_struct_array->GetFieldByName("trackStates"));    auto arr_trackStates_val = std::static_pointer_cast<arrow::StructArray>(arr_trackStates->values());    auto arr_trackStates_val_location = std::static_pointer_cast<arrow::Int32Array>(arr_trackStates_val->GetFieldByName("location"));
    auto arr_trackStates_val_D0 = std::static_pointer_cast<arrow::FloatArray>(arr_trackStates_val->GetFieldByName("D0"));
    auto arr_trackStates_val_phi = std::static_pointer_cast<arrow::FloatArray>(arr_trackStates_val->GetFieldByName("phi"));
    auto arr_trackStates_val_omega = std::static_pointer_cast<arrow::FloatArray>(arr_trackStates_val->GetFieldByName("omega"));
    auto arr_trackStates_val_Z0 = std::static_pointer_cast<arrow::FloatArray>(arr_trackStates_val->GetFieldByName("Z0"));
    auto arr_trackStates_val_tanLambda = std::static_pointer_cast<arrow::FloatArray>(arr_trackStates_val->GetFieldByName("tanLambda"));
    auto arr_trackStates_val_time = std::static_pointer_cast<arrow::FloatArray>(arr_trackStates_val->GetFieldByName("time"));
    auto arr_trackStates_val_referencePoint = std::static_pointer_cast<arrow::StructArray>(arr_trackStates_val->GetFieldByName("referencePoint"));
    auto arr_trackStates_val_referencePoint_x = std::static_pointer_cast<arrow::FloatArray>(arr_trackStates_val_referencePoint->GetFieldByName("x"));

    auto arr_trackStates_val_referencePoint_y = std::static_pointer_cast<arrow::FloatArray>(arr_trackStates_val_referencePoint->GetFieldByName("y"));

    auto arr_trackStates_val_referencePoint_z = std::static_pointer_cast<arrow::FloatArray>(arr_trackStates_val_referencePoint->GetFieldByName("z"));

    auto arr_trackStates_val_covMatrix = std::static_pointer_cast<arrow::StructArray>(arr_trackStates_val->GetFieldByName("covMatrix"));
    auto arr_trackStates_val_covMatrix_values = std::static_pointer_cast<arrow::FixedSizeListArray>(arr_trackStates_val_covMatrix->GetFieldByName("values"));
    auto arr_trackStates_val_covMatrix_values_val = std::static_pointer_cast<arrow::FloatArray>(arr_trackStates_val_covMatrix_values->values());

    auto* vecMemberBuffer_trackStates = podio::CollectionReadBuffers::asVector<edm4hep::TrackState>(buffers->vectorMembers->at(1).second);
    vecMemberBuffer_trackStates->reserve(arr_trackStates->values()->length());    auto arr_dxQuantities = std::static_pointer_cast<arrow::ListArray>(raw_data_struct_array->GetFieldByName("dxQuantities"));    auto arr_dxQuantities_val = std::static_pointer_cast<arrow::StructArray>(arr_dxQuantities->values());    auto arr_dxQuantities_val_type = std::static_pointer_cast<arrow::Int32Array>(arr_dxQuantities_val->GetFieldByName("type"));
    auto arr_dxQuantities_val_value = std::static_pointer_cast<arrow::FloatArray>(arr_dxQuantities_val->GetFieldByName("value"));
    auto arr_dxQuantities_val_error = std::static_pointer_cast<arrow::FloatArray>(arr_dxQuantities_val->GetFieldByName("error"));
    auto* vecMemberBuffer_dxQuantities = podio::CollectionReadBuffers::asVector<edm4hep::Quantity>(buffers->vectorMembers->at(2).second);
    vecMemberBuffer_dxQuantities->reserve(arr_dxQuantities->values()->length());    auto rel_arr_trackerHits = std::static_pointer_cast<arrow::ListArray>(raw_data_struct_array->GetFieldByName("trackerHits"));
    auto rel_struct_arr_trackerHits = std::static_pointer_cast<arrow::StructArray>(rel_arr_trackerHits->values());
    auto rel_collId_arr_trackerHits = std::static_pointer_cast<arrow::UInt32Array>(rel_struct_arr_trackerHits->field(0));
    auto rel_index_arr_trackerHits = std::static_pointer_cast<arrow::Int32Array>(rel_struct_arr_trackerHits->field(1));
    auto* refVec_trackerHits = buffers->references->at(0).get();
    refVec_trackerHits->reserve(rel_arr_trackerHits->values()->length());    auto rel_arr_tracks = std::static_pointer_cast<arrow::ListArray>(raw_data_struct_array->GetFieldByName("tracks"));
    auto rel_struct_arr_tracks = std::static_pointer_cast<arrow::StructArray>(rel_arr_tracks->values());
    auto rel_collId_arr_tracks = std::static_pointer_cast<arrow::UInt32Array>(rel_struct_arr_tracks->field(0));
    auto rel_index_arr_tracks = std::static_pointer_cast<arrow::Int32Array>(rel_struct_arr_tracks->field(1));
    auto* refVec_tracks = buffers->references->at(1).get();
    refVec_tracks->reserve(rel_arr_tracks->values()->length());
    for (size_t i = 0; i < collection_size; ++i) {      (*dataVec)[i].type = arr_type->Value(i);
      (*dataVec)[i].chi2 = arr_chi2->Value(i);
      (*dataVec)[i].ndf = arr_ndf->Value(i);
      (*dataVec)[i].dEdx = arr_dEdx->Value(i);
      (*dataVec)[i].dEdxError = arr_dEdxError->Value(i);
      (*dataVec)[i].radiusOfInnermostHit = arr_radiusOfInnermostHit->Value(i);
      {
        size_t count_begin = vecMemberBuffer_subdetectorHitNumbers->size();
        // Get start and end offsets of the vector/list element in the flat values array
        int32_t start_offset = arr_subdetectorHitNumbers->value_offset(i);
        int32_t end_offset = arr_subdetectorHitNumbers->value_offset(i + 1);
        for (int32_t j = start_offset; j < end_offset; ++j) {
vecMemberBuffer_subdetectorHitNumbers->push_back(arr_subdetectorHitNumbers_val->Value(j));
}
        (*dataVec)[i].subdetectorHitNumbers_begin = count_begin;
        (*dataVec)[i].subdetectorHitNumbers_end = vecMemberBuffer_subdetectorHitNumbers->size();
      }      {
        size_t count_begin = vecMemberBuffer_trackStates->size();
        // Get start and end offsets of the vector/list element in the flat values array
        int32_t start_offset = arr_trackStates->value_offset(i);
        int32_t end_offset = arr_trackStates->value_offset(i + 1);
        for (int32_t j = start_offset; j < end_offset; ++j) {
edm4hep::TrackState val;
val.location = arr_trackStates_val_location->Value(j);

val.D0 = arr_trackStates_val_D0->Value(j);

val.phi = arr_trackStates_val_phi->Value(j);

val.omega = arr_trackStates_val_omega->Value(j);

val.Z0 = arr_trackStates_val_Z0->Value(j);

val.tanLambda = arr_trackStates_val_tanLambda->Value(j);

val.time = arr_trackStates_val_time->Value(j);

val.referencePoint.x = arr_trackStates_val_referencePoint_x->Value(j);

          val.referencePoint.y = arr_trackStates_val_referencePoint_y->Value(j);

          val.referencePoint.z = arr_trackStates_val_referencePoint_z->Value(j);


for (size_t arr_i = 0; arr_i < 21; ++arr_i) {
          val.covMatrix.values[arr_i] = arr_trackStates_val_covMatrix_values_val->Value((arr_trackStates_val_covMatrix_values->offset() + j) * 21 + arr_i);
          }


vecMemberBuffer_trackStates->push_back(val);
}
        (*dataVec)[i].trackStates_begin = count_begin;
        (*dataVec)[i].trackStates_end = vecMemberBuffer_trackStates->size();
      }      {
        size_t count_begin = vecMemberBuffer_dxQuantities->size();
        // Get start and end offsets of the vector/list element in the flat values array
        int32_t start_offset = arr_dxQuantities->value_offset(i);
        int32_t end_offset = arr_dxQuantities->value_offset(i + 1);
        for (int32_t j = start_offset; j < end_offset; ++j) {
edm4hep::Quantity val;
val.type = arr_dxQuantities_val_type->Value(j);

val.value = arr_dxQuantities_val_value->Value(j);

val.error = arr_dxQuantities_val_error->Value(j);

vecMemberBuffer_dxQuantities->push_back(val);
}
        (*dataVec)[i].dxQuantities_begin = count_begin;
        (*dataVec)[i].dxQuantities_end = vecMemberBuffer_dxQuantities->size();
      }      {
        size_t count_begin = refVec_trackerHits->size();
        // Get start and end offsets of the relation list in the flat values array
        int32_t start_offset = rel_arr_trackerHits->value_offset(i);
        int32_t end_offset = rel_arr_trackerHits->value_offset(i + 1);
        for (int32_t j = start_offset; j < end_offset; ++j) {
          refVec_trackerHits->push_back(podio::ObjectID{rel_index_arr_trackerHits->Value(j), rel_collId_arr_trackerHits->Value(j)});
        }
        (*dataVec)[i].trackerHits_begin = count_begin;
        (*dataVec)[i].trackerHits_end = refVec_trackerHits->size();
      }      {
        size_t count_begin = refVec_tracks->size();
        // Get start and end offsets of the relation list in the flat values array
        int32_t start_offset = rel_arr_tracks->value_offset(i);
        int32_t end_offset = rel_arr_tracks->value_offset(i + 1);
        for (int32_t j = start_offset; j < end_offset; ++j) {
          refVec_tracks->push_back(podio::ObjectID{rel_index_arr_tracks->Value(j), rel_collId_arr_tracks->Value(j)});
        }
        (*dataVec)[i].tracks_begin = count_begin;
        (*dataVec)[i].tracks_end = refVec_tracks->size();
      }    }
    return buffers;
  });
  convRegistry.registerReader("extension::SenseWireHit", [](std::shared_ptr<arrow::Array> array, int64_t rowIndex, bool isSubset, podio::SchemaVersionT version) -> std::optional<podio::CollectionReadBuffers> {
    auto buffers = podio::CollectionBufferFactory::instance().createBuffers("extension::SenseWireHitCollection", version, isSubset);
    if (!buffers) {
      return std::nullopt;
    }

    auto collection_frames_list = std::static_pointer_cast<arrow::ListArray>(array);
    auto frame_collection_slice = collection_frames_list->value_slice(rowIndex);
    auto raw_data_struct_array = std::static_pointer_cast<arrow::StructArray>(frame_collection_slice);

    if (!buffers->data) {
      auto collId_arr = std::static_pointer_cast<arrow::UInt32Array>(raw_data_struct_array->field(0));
      auto index_arr = std::static_pointer_cast<arrow::Int32Array>(raw_data_struct_array->field(1));
      auto* refVec = buffers->references->at(0).get();
      size_t collection_size = raw_data_struct_array->length();
      refVec->reserve(collection_size);
      for (size_t i = 0; i < collection_size; ++i) {
        refVec->push_back(podio::ObjectID{index_arr->Value(i), collId_arr->Value(i)});
      }
      return buffers;
    }

    auto* dataVec = buffers->dataAsVector<extension::SenseWireHitData>();
    size_t collection_size = raw_data_struct_array->length();
    dataVec->resize(collection_size);    auto arr_cellID = std::static_pointer_cast<arrow::UInt64Array>(raw_data_struct_array->GetFieldByName("cellID"));
    auto arr_type = std::static_pointer_cast<arrow::Int32Array>(raw_data_struct_array->GetFieldByName("type"));
    auto arr_quality = std::static_pointer_cast<arrow::Int32Array>(raw_data_struct_array->GetFieldByName("quality"));
    auto arr_time = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("time"));
    auto arr_eDep = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("eDep"));
    auto arr_eDepError = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("eDepError"));
    auto arr_wireStereoAngle = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("wireStereoAngle"));
    auto arr_wireAzimuthalAngle = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("wireAzimuthalAngle"));
    auto arr_position = std::static_pointer_cast<arrow::StructArray>(raw_data_struct_array->GetFieldByName("position"));
    auto arr_position_x = std::static_pointer_cast<arrow::DoubleArray>(arr_position->GetFieldByName("x"));

    auto arr_position_y = std::static_pointer_cast<arrow::DoubleArray>(arr_position->GetFieldByName("y"));

    auto arr_position_z = std::static_pointer_cast<arrow::DoubleArray>(arr_position->GetFieldByName("z"));

    auto arr_positionAlongWireError = std::static_pointer_cast<arrow::DoubleArray>(raw_data_struct_array->GetFieldByName("positionAlongWireError"));
    auto arr_distanceToWire = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("distanceToWire"));
    auto arr_distanceToWireError = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("distanceToWireError"));
    auto arr_nElectrons = std::static_pointer_cast<arrow::ListArray>(raw_data_struct_array->GetFieldByName("nElectrons"));    auto arr_nElectrons_val = std::static_pointer_cast<arrow::UInt16Array>(arr_nElectrons->values());    auto* vecMemberBuffer_nElectrons = podio::CollectionReadBuffers::asVector<std::uint16_t>(buffers->vectorMembers->at(0).second);
    vecMemberBuffer_nElectrons->reserve(arr_nElectrons->values()->length());
    for (size_t i = 0; i < collection_size; ++i) {      (*dataVec)[i].cellID = arr_cellID->Value(i);
      (*dataVec)[i].type = arr_type->Value(i);
      (*dataVec)[i].quality = arr_quality->Value(i);
      (*dataVec)[i].time = arr_time->Value(i);
      (*dataVec)[i].eDep = arr_eDep->Value(i);
      (*dataVec)[i].eDepError = arr_eDepError->Value(i);
      (*dataVec)[i].wireStereoAngle = arr_wireStereoAngle->Value(i);
      (*dataVec)[i].wireAzimuthalAngle = arr_wireAzimuthalAngle->Value(i);
      (*dataVec)[i].position.x = arr_position_x->Value(i);

      (*dataVec)[i].position.y = arr_position_y->Value(i);

      (*dataVec)[i].position.z = arr_position_z->Value(i);

      (*dataVec)[i].positionAlongWireError = arr_positionAlongWireError->Value(i);
      (*dataVec)[i].distanceToWire = arr_distanceToWire->Value(i);
      (*dataVec)[i].distanceToWireError = arr_distanceToWireError->Value(i);
      {
        size_t count_begin = vecMemberBuffer_nElectrons->size();
        // Get start and end offsets of the vector/list element in the flat values array
        int32_t start_offset = arr_nElectrons->value_offset(i);
        int32_t end_offset = arr_nElectrons->value_offset(i + 1);
        for (int32_t j = start_offset; j < end_offset; ++j) {
vecMemberBuffer_nElectrons->push_back(arr_nElectrons_val->Value(j));
}
        (*dataVec)[i].nElectrons_begin = count_begin;
        (*dataVec)[i].nElectrons_end = vecMemberBuffer_nElectrons->size();
      }    }
    return buffers;
  });
  convRegistry.registerReader("extension::SenseWireHitSimTrackerHitLink", [](std::shared_ptr<arrow::Array> array, int64_t rowIndex, bool isSubset, podio::SchemaVersionT version) -> std::optional<podio::CollectionReadBuffers> {
    auto buffers = podio::CollectionBufferFactory::instance().createBuffers("extension::SenseWireHitSimTrackerHitLinkCollection", version, isSubset);
    if (!buffers) {
      return std::nullopt;
    }

    auto collection_frames_list = std::static_pointer_cast<arrow::ListArray>(array);
    auto frame_collection_slice = collection_frames_list->value_slice(rowIndex);
    auto raw_data_struct_array = std::static_pointer_cast<arrow::StructArray>(frame_collection_slice);

    if (!buffers->data) {
      auto collId_arr = std::static_pointer_cast<arrow::UInt32Array>(raw_data_struct_array->field(0));
      auto index_arr = std::static_pointer_cast<arrow::Int32Array>(raw_data_struct_array->field(1));
      auto* refVec = buffers->references->at(0).get();
      size_t collection_size = raw_data_struct_array->length();
      refVec->reserve(collection_size);
      for (size_t i = 0; i < collection_size; ++i) {
        refVec->push_back(podio::ObjectID{index_arr->Value(i), collId_arr->Value(i)});
      }
      return buffers;
    }

    auto* dataVec = buffers->dataAsVector<extension::SenseWireHitSimTrackerHitLinkData>();
    size_t collection_size = raw_data_struct_array->length();
    dataVec->resize(collection_size);    auto arr_weight = std::static_pointer_cast<arrow::FloatArray>(raw_data_struct_array->GetFieldByName("weight"));
    auto rel_struct_arr_from = std::static_pointer_cast<arrow::StructArray>(raw_data_struct_array->GetFieldByName("from"));
    auto rel_collId_arr_from = std::static_pointer_cast<arrow::UInt32Array>(rel_struct_arr_from->field(0));
    auto rel_index_arr_from = std::static_pointer_cast<arrow::Int32Array>(rel_struct_arr_from->field(1));
    auto* refVec_from = buffers->references->at(0).get();
    refVec_from->reserve(collection_size);    auto rel_struct_arr_to = std::static_pointer_cast<arrow::StructArray>(raw_data_struct_array->GetFieldByName("to"));
    auto rel_collId_arr_to = std::static_pointer_cast<arrow::UInt32Array>(rel_struct_arr_to->field(0));
    auto rel_index_arr_to = std::static_pointer_cast<arrow::Int32Array>(rel_struct_arr_to->field(1));
    auto* refVec_to = buffers->references->at(1).get();
    refVec_to->reserve(collection_size);
    for (size_t i = 0; i < collection_size; ++i) {      (*dataVec)[i].weight = arr_weight->Value(i);
      {
        refVec_from->push_back(podio::ObjectID{rel_index_arr_from->Value(i), rel_collId_arr_from->Value(i)});
      }      {
        refVec_to->push_back(podio::ObjectID{rel_index_arr_to->Value(i), rel_collId_arr_to->Value(i)});
      }    }
    return buffers;
  });
  return true;
}

static const bool readersRegistered = registerBufferReaders();

} // namespace

} // namespace extension::arrow_io
