// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_extension_H
#define EXTENSION_extension_H

#include "podio/utilities/TypeHelpers.h"

#include "extension/DriftChamberDigiLocalCollection.h"
#include "extension/DriftChamberDigiCollection.h"
#include "extension/MCRecoDriftChamberDigiAssociationCollection.h"
#include "extension/DriftChamberDigiV2Collection.h"
#include "extension/MCRecoDriftChamberDigiV2AssociationCollection.h"
#include "extension/TrackCollection.h"
#include "extension/SenseWireHitCollection.h"
#include "extension/SenseWireHitSimTrackerHitLinkCollection.h"

namespace extension {

using extensionDataTypes = podio::utils::TypeList<
  extension::DriftChamberDigiLocal, extension::DriftChamberDigi, extension::MCRecoDriftChamberDigiAssociation, extension::DriftChamberDigiV2, extension::MCRecoDriftChamberDigiV2Association, extension::Track, extension::SenseWireHit, extension::SenseWireHitSimTrackerHitLink
>;

using extensionDataCollectionTypes = podio::utils::TypeList<
  extension::DriftChamberDigiLocalCollection, extension::DriftChamberDigiCollection, extension::MCRecoDriftChamberDigiAssociationCollection, extension::DriftChamberDigiV2Collection, extension::MCRecoDriftChamberDigiV2AssociationCollection, extension::TrackCollection, extension::SenseWireHitCollection, extension::SenseWireHitSimTrackerHitLinkCollection
>;

using extensionLinkTypes = podio::utils::TypeList<
  
>;

using extensionLinkCollectionTypes = podio::utils::TypeList<
>;

using extensionInterfaceTypes = podio::utils::TypeList<
  extension::TrackerHit, extension::TrackerHit_dev
>;

}

#endif // EXTENSION_extension_H
