#!/usr/bin/env python
# MAP from java change map if is needed
import json,yaml,ast

maptext ='''UNDEFINED ( 0, "UNDEF"),
    BMT       ( 1, "BMT"),    
    BST       ( 2, "BST"),
    CND       ( 3, "CND"),
    CTOF      ( 4, "CTOF"),
    CVT       ( 5, "CVT"),
    DC        ( 6, "DC"),
    ECAL      ( 7, "ECAL"),
    FMT       ( 8, "FMT"),
    FT        ( 9, "FT"),
    FTCAL     (10, "FTCAL"),
    FTHODO    (11, "FTHODO"),
    FTOF      (12, "FTOF"),
    FTTRK     (13, "FTTRK"),
    HTCC      (15, "HTCC"),
    LTCC      (16, "LTCC"),
    RF        (17, "RF"),
    RICH      (18, "RICH"),
    RTPC      (19, "RTPC"),
    HEL       (20, "HEL"),
    BAND      (21, "BAND"),
    ECIN      (110, "ECIN"),
    ECOUT     (111, "ECOUT"),
    ECTOT     (112, "ECTOT"),
    LAC       (113, "LAC"),
    SC        (114, "SC"),
    CC        (115, "CC")
'''
maptext = (map(lambda x: x.strip()[10:].replace(",",":").replace("(","").replace(")","") ,maptext.strip().split('),')))
maptext = reduce(lambda x,y: x + "," + y,maptext)
maptext = "{" + maptext + "}"
#d =json.loads(maptext)
d=ast.literal_eval(maptext)
for (k,v) in zip(d.keys(),d.values()):
    print 'detectorType["' + v + '"] = '+str(k) + ";"
    
 
