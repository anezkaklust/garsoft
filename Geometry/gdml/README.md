Description of the gdml files

MPD Volumes by descending order are:

- volWorld
---- rockBox_lv (if in Hall)
---- volDetEnclosure
------- volMPD
--------- volNDHPgTPC
------------ volGArTPC
--------------- volTPCChamber
------------------ volTPCGas
------------ volPVBarrel
------------ volPVEndcap
------------ volBarrelECal
--------------- BarrelECal_stave%i_module%i_vol (staves 1 to 8 and modules 1 to 5)
------------ volEndcapECal
--------------- EndcapECal_stave%i_module%i_vol (staves 1 to 4 and modules 0 and 6)
--------- volMagnet
------------ volMagnet_Coil%i (if Coils)
--------- volYokeBarrel (if Yoke)
------------ YokeBarrel_stave%i_module%i (staves 1 to 8 and modules 1 and 2)
--------- volYokeEndcap (if Yoke)

=============================================================

MPD only in world volume (located in folder MPD_Standalone)

- ND_Strawman_Concept_v09.gdml

Legacy MPD concept at the ND CDR-lite time

- MPD_ECalOctagon_60l_UniformMagnet.gdml

Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and a uniform cylinder magnet of Aluminium of around 100t.

- MPD_ECalDodecagon_80l_UniformMagnet.gdml

MPD with HPgTPC and PV, an ECAL in Dodecagonal shape (12 sides) with 80 layers in Barrel and 60 layers in the Endcap of 2 mm Cu, 5 mm Sc and a uniform cylinder magnet of Aluminium of around 100t.

- MPD_ECalOctagon_60l_2Coils.gdml

Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and 2 Coils Helmholtz magnet

- MPD_ECalOctagon_60l_4Coils.gdml

Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and 4 Coils Helmholtz magnet

- MPD_ECalOctagon_60l_5Coils.gdml

Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and 5 Coils Helmholtz magnet

- MPD_ECalOctagon_60l_5Coils_ThinUpstreamECAL.gdml

Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and 5 Coils Helmholtz magnet. The upstream part of the ECAL is thinner and consists of 8 HG layers

- MPD_ECalOctagon_60l_5Coils_noUpstreamECAL.gdml

Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and 5 Coils Helmholtz magnet. No upstream ECAL

- MPD_ECalOctagon_60l_SPY.gdml (https://indico.fnal.gov/event/21340/session/5/contribution/40/material/slides/0.pdf)

Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and SPY magnet/PRY design.

- MPD_ECalOctagon_60l_SPY_ThinUpstreamECAL.gdml

Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and SPY magnet/PRY design. The upstream part of the ECAL is thinner and consists of 8 HG layers

- MPD_ECalOctagon_60l_SPY_noUpstreamECAL.gdml

Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and SPY magnet/PRY design. No Upstream ECAL

- MPD_ECalOctagon_60l_SPY_wMuID.gdml

Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and SPY magnet/PRY design instrumented with 3 layers for muID (5 cm Iron + 1.67 cm Sc).

- MPD_ECalOctagon_60l_SPY_wMuID_ThinUpstreamECAL.gdml

Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and SPY magnet/PRY design instrumented with 3 layers for muID (5 cm Iron + 1.67 cm Sc). The upstream part of the ECAL is thinner and consists of 8 HG layers

- MPD_ECalOctagon_60l_SPY_wMuID_noUpstreamECAL.gdml

Baseline MPD with HPgTPC and PV, an ECAL in Octagonal shape (8 sides) with 60 layers in Barrel and Endcap of 2 mm Cu, 5 mm Sc and SPY magnet/PRY design instrumented with 3 layers for muID (5 cm Iron + 1.67 cm Sc). No Upstream ECAL

=============================================================

MPD only in ND Hall (located in folder nd_hall_mpd_only)

These geometries are the same as above (MPD standalone) however the MPD is placed in the ND Hall

- nd_hall_mpd_only_ECalDodecagon_80l_UniformMagnet.gdml

- nd_hall_mpd_only_ECalOctagon_60l_2Coils.gdml

- nd_hall_mpd_only_ECalOctagon_60l_4Coils.gdml

- nd_hall_mpd_only_ECalOctagon_60l_5Coils.gdml

- nd_hall_mpd_only_ECalOctagon_60l_5Coils_noUpstreamECAL.gdml

- nd_hall_mpd_only_ECalOctagon_60l_5Coils_ThinUpstreamECAL.gdml

- nd_hall_mpd_only_ECalOctagon_60l_SPY.gdml

- nd_hall_mpd_only_ECalOctagon_60l_SPY_noUpstreamECAL.gdml

- nd_hall_mpd_only_ECalOctagon_60l_SPY_ThinUpstreamECAL.gdml

- nd_hall_mpd_only_ECalOctagon_60l_SPY_wMuID.gdml

- nd_hall_mpd_only_ECalOctagon_60l_SPY_wMuID_ThinUpstreamECAL.gdml

- nd_hall_mpd_only_ECalOctagon_60l_UniformMagnet.gdml

- nd_hall_mpd_only_ECalOctagon_60l_SPY_wMuID_noUpstreamECAL.gdml

=============================================================

MPD with LAr in ND Hall (located in folder nd_hall_mpd_lar_only)

These geometries are the same as above (MPD standalone) however the MPD is placed in the ND Hall along side the LArTPC

- nd_hall_mpd_lar_ECalOctagon_60l_UniformMagnet.gdml

- nd_hall_mpd_lar_ECalOctagon_60l_SPY_wMuID.gdml

- nd_hall_mpd_lar_ECalOctagon_60l_SPY.gdml

- nd_hall_mpd_lar_ECalOctagon_60l_5Coils_ThinUpstreamECAL.gdml

- nd_hall_mpd_lar_ECalOctagon_60l_5Coils_noUpstreamECAL.gdml

- nd_hall_mpd_lar_ECalOctagon_60l_5Coils.gdml

- nd_hall_mpd_lar_ECalOctagon_60l_SPY_wMuID_ThinUpstreamECAL.gdml

- nd_hall_mpd_lar_ECalOctagon_60l_SPY_ThinUpstreamECAL.gdml

- nd_hall_mpd_lar_ECalOctagon_60l_SPY_noUpstreamECAL.gdml

- nd_hall_mpd_lar_ECalOctagon_60l_SPY_wMuID_noUpstreamECAL.gdml

=============================================================

Full ND Complex (located in folder nd_hall_all_dets)

These geometries are the same as above (MPD standalone) however the MPD is placed in the ND Hall along side the LArTPC and SAND

- nd_hall_all_dets_ECalOctagon_60l_UniformMagnet.gdml

- nd_hall_all_dets_ECalOctagon_60l_SPY.gdml

- nd_hall_all_dets_ECalOctagon_60l_SPY_wMuID_ThinUpstreamECAL.gdml

- nd_hall_all_dets_ECalOctagon_60l_SPY_wMuID.gdml

- nd_hall_all_dets_ECalOctagon_60l_SPY_ThinUpstreamECAL.gdml

- nd_hall_all_dets_ECalOctagon_60l_SPY_noUpstreamECAL.gdml

- nd_hall_all_dets_ECalOctagon_60l_5Coils_ThinUpstreamECAL.gdml

- nd_hall_all_dets_ECalOctagon_60l_5Coils_noUpstreamECAL.gdml

- nd_hall_all_dets_ECalOctagon_60l_5Coils.gdml

- nd_hall_all_dets_ECalOctagon_60l_SPY_wMuID_noUpstreamECAL.gdml

=============================================================
