import numpy as np
import json
from scipy.signal import argrelextrema
from pyaedt import Hfss, Desktop

class HPTHFSS(object):
    def __init__(self, kwargs, ansys_handles, index=None):
        super(HPTHFSS, self).__init__()
        self.kwargs = kwargs	# 参数
        self.index = index		# 第几次仿真
        self.obj_names = {"RegularPolyhedron":[{"Board": []}, 1],
                          "Polygon":[{"Coil_out":[], "Coil_in":[]},1],
                          "Rectangle":[{"Gap":[], "Lead":[]},1],
                          "Polyline": [{"Connector": []}, 1],
                        #   "Exitation":[{"Exitation":[]},2],
                          "SecondPCB":[{"Second":[]},1]
                          } # EXTENSION
        self.subtracted_names = []
        self.boundary = []
        self.PCBspace = 10 # mm
        self.frequencyrange = [1, 5]  # 粗略扫频的频率范围, 要看看SRF在什么地方，如果比较小的话感觉不太需要扫描两次
        self.res_path = kwargs["res_path"]
        self.report_files = []	# 保存文件夹名字
        self.ansys_handles = ansys_handles # ansys_handles := (oDesktop, oProject, oDesign, oEditor)
        self.compute_rad()
    
    def compute_rad(self):  
        wins = self.kwargs["wins"]
        wouts = self.kwargs["wouts"]
        
        # 初始化外圈内径和外径
        self.radins1 = [self.kwargs["rout_out"] - wouts[0]]     # 初始化外圈内径
        self.radouts1 = [self.kwargs["rout_out"]]               # 初始化外圈外径        

        # 遍历 wouts 列表，跳过第一个元素
        for i in range(1, len(wouts)):
            current_w = wouts[i]          # 当前层的厚度
            previous_w = wouts[i - 1]     # 上一层的厚度

            # 计算新的内径和外径
            radin = self.radins1[-1] - self.kwargs["space"] - current_w
            radout = self.radouts1[-1] - self.kwargs["space"] - previous_w

            # 添加到列表中
            self.radins1.append(radin)
            self.radouts1.append(radout)
            
        # 初始化内圈内径和外径
        self.radins2 = [self.kwargs["rout_in"] - wins[0]]   # 初始化内圈内径
        self.radouts2 = [self.kwargs["rout_in"]]  # 初始化内圈外径
        
        # 遍历 wins 列表，跳过第一个元素
        for i in range(1, len(wins)):
            current_w = wins[i]          # 当前层的厚度
            previous_w = wins[i - 1]     # 上一层的厚度

            # 计算新的内径和外径
            radin = self.radins2[-1] - self.kwargs["space"] - current_w
            radout = self.radouts2[-1] - self.kwargs["space"] - previous_w

            # 添加到列表中
            self.radins2.append(radin)
            self.radouts2.append(radout)

    # 在材料库中创建材料PTFE
    def create_PTFE(self):
        oProject = self.ansys_handles[-3]
        oDefinitionManager = oProject.GetDefinitionManager()
        oDefinitionManager.AddMaterial(
            [
                "NAME:PTFE",
                "CoordinateSystemType:=", "Cartesian",
                "BulkOrSurfaceType:="	, 1,
                [
                    "NAME:PhysicsTypes",
                    "set:="			, ["Electromagnetic","Thermal","Structural"]
                ],
                [
                    "NAME:AttachedData",
                    [
                        "NAME:MatAppearanceData",
                        "property_data:="	, "appearance_data",
                        "Red:="			, 27,
                        "Green:="		, 110,
                        "Blue:="		, 76
                    ]
                ],
                "permittivity:="	, "2.94",
                "dielectric_loss_tangent:=", "0.00016",
                "thermal_conductivity:=", "0.294",
                "mass_density:="	, "2200",
                "specific_heat:="	, "1150",
                "youngs_modulus:="	, "11000000000",
                "poissons_ratio:="	, "0.28",
                "thermal_expansion_coefficient:=", "1.5e-05"
            ])
        
    # 绘制PCB基板，PCB尺寸固定
    def create_PCB(self):   
        oEditor = self.ansys_handles[-1]
        
        RegularPolyhedron_name = oEditor.CreateRegularPolyhedron(
            [
                "NAME:PolyhedronParameters",
                "XCenter:="		, "0mm",
                "YCenter:="		, "0mm",
                "ZCenter:="		, "0mm",
                "XStart:="		, "0mm",
                "YStart:="		, "160mm",	# EXTENSION
                "ZStart:="		, "0mm",
                "Height:="		, "-1.52mm",
                "NumSides:="		, "58",
                "WhichAxis:="		, "Z"
            ], 
            [
                "NAME:Attributes",
                "Name:="		, "RegularPolyhedron1",
                "Flags:="		, "",
                "Color:="		, "(143 175 143)",
                "Transparency:="	, 0,
                "PartCoordinateSystem:=", "Global",
                "UDMId:="		, "",
                "MaterialValue:="	, "\"PTFE\"",
                "SurfaceMaterialValue:=", "\"\"",
                "SolveInside:="		, True,
                "ShellElement:="	, False,
                "ShellElementThickness:=", "0mm",
                "ReferenceTemperature:=", "20cel",
                "IsMaterialEditable:="	, True,
                "UseMaterialAppearance:=", False,
                "IsLightweight:="	, False
            ])
        self.obj_names["RegularPolyhedron"][0]["Board"].append(RegularPolyhedron_name)
        self.obj_names["RegularPolyhedron"][1] += 1
        
    def create_coils(self):
        radins1, radouts1, radins2, radouts2 = self.radins1, self.radouts1, self.radins2, self.radouts2
        n1, n2 = self.kwargs["n1"], self.kwargs["n2"]
        oEditor = self.ansys_handles[-1]
        # 先画外圈
        for i in range(1, 2*n1+1):
            polygon_name = f"Polygon{i}"
            index = (i - 1) // 2  # 计算当前的索引，每两次循环后索引增加 1
            if i % 2 == 1:  # 如果 i 是奇数
               radius = radouts1[index]
            else:  # 如果 i 是偶数
                radius = radins1[index]
            oEditor.CreateRegularPolygon(  
                [
                    "NAME:RegularPolygonParameters",
                    "IsCovered:="		, True,
                    "XCenter:="		, "0mm",
                    "YCenter:="		, "0mm",
                    "ZCenter:="		, "0mm",
                    "XStart:="		, "0mm",
                    "YStart:="		, f"{radius}mm",  #EXTENSION
                    "ZStart:="		, "0mm",
                    "NumSides:="		, "58",
                    "WhichAxis:="		, "Z"
                ], 
                [
                    "NAME:Attributes",
                    "Name:="		, polygon_name,   
                    "Flags:="		, "",
                    "Color:="		, "(143 175 143)",
                    "Transparency:="	, 0,
                    "PartCoordinateSystem:=", "Global",
                    "UDMId:="		, "",
                    "MaterialValue:="	, "\"vacuum\"",
                        "SurfaceMaterialValue:=", "\"\"",
                    "SolveInside:="		, True,
                    "ShellElement:="	, False,
                    "ShellElementThickness:=", "0mm",
                    "ReferenceTemperature:=", "20cel",
                    "IsMaterialEditable:="	, True,
                    "UseMaterialAppearance:=", False,
                    "IsLightweight:="	, False
                ])        
            self.obj_names["Polygon"][0]["Coil_out"].append(polygon_name)
            self.obj_names["Polygon"][1] += 1   # 每新建一个多边形，这个值就+1，创建了Polygon1对应这个值是2
            
            if i % 2 == 0:  # 如果i是偶数，就Polygon{i-1}-Polygon{i}
                oEditor.Subtract(
                    [
                        "NAME:Selections",
                        "Blank Parts:="		, f"Polygon{self.obj_names['Polygon'][1]-2}",
                        "Tool Parts:="		, f"Polygon{self.obj_names['Polygon'][1]-1}"
                    ], 
                    [
                        "NAME:SubtractParameters",
                        "KeepOriginals:="	, False,
                        "TurnOnNBodyBoolean:="	, True
                    ])
                self.subtracted_names.append(f"Polygon{self.obj_names['Polygon'][1]-1}")    # 记录被减去的部分的名字

        for i in range(1, 2*n2+1):
            polygon_name = f"Polygon{2*n1+i}"
            index = (i-1)//2
            if i % 2 == 1:
                radius = radouts2[index]
            else:
                radius = radins2[index]
            oEditor.CreateRegularPolygon(  
                [
                    "NAME:RegularPolygonParameters",
                    "IsCovered:="		, True,
                    "XCenter:="		, "0mm",
                    "YCenter:="		, "0mm",
                    "ZCenter:="		, "0mm",
                    "XStart:="		, "0mm",
                    "YStart:="		, f"{radius}mm",  #EXTENSION
                    "ZStart:="		, "0mm",
                    "NumSides:="		, "58",
                    "WhichAxis:="		, "Z"
                ], 
                [
                    "NAME:Attributes",
                    "Name:="		, polygon_name,   
                    "Flags:="		, "",
                    "Color:="		, "(143 175 143)",
                    "Transparency:="	, 0,
                    "PartCoordinateSystem:=", "Global",
                    "UDMId:="		, "",
                    "MaterialValue:="	, "\"vacuum\"",
                        "SurfaceMaterialValue:=", "\"\"",
                    "SolveInside:="		, True,
                    "ShellElement:="	, False,
                    "ShellElementThickness:=", "0mm",
                    "ReferenceTemperature:=", "20cel",
                    "IsMaterialEditable:="	, True,
                    "UseMaterialAppearance:=", False,
                    "IsLightweight:="	, False
                ])
            self.obj_names["Polygon"][0]["Coil_in"].append(polygon_name)
            self.obj_names["Polygon"][1] += 1            

            if i % 2 == 0:  # 如果i是偶数，就Polygon{i-1}-Polygon{i}
                oEditor.Subtract(
                    [
                        "NAME:Selections",
                        "Blank Parts:="		, f"Polygon{self.obj_names['Polygon'][1]-2}",
                        "Tool Parts:="		, f"Polygon{self.obj_names['Polygon'][1]-1}"
                    ], 
                    [
                        "NAME:SubtractParameters",
                        "KeepOriginals:="	, False,
                        "TurnOnNBodyBoolean:="	, True
                    ])
                self.subtracted_names.append(f"Polygon{self.obj_names['Polygon'][1]-1}")
            
    def create_gap(self):
        rout_out, oEditor = self.kwargs["rout_out"], self.ansys_handles[-1]
        rectangle_name = oEditor.CreateRectangle(
            [
            "NAME:RectangleParameters",
            "IsCovered:="		, True,
            "XStart:="		, "0mm", 
            "YStart:="		, "-4mm",   # fixed
            "ZStart:="		, "0mm",
            "Width:="		, f"-{rout_out+10}mm",  #EXTENSION
            "Height:="		, "8mm",    # fixed
            "WhichAxis:="		, "Z"
            ], 
            [
            "NAME:Attributes",
            "Name:="		, "Rectangle1",
            "Flags:="		, "",
            "Color:="		, "(143 175 143)",
            "Transparency:="	, 0,
            "PartCoordinateSystem:=", "Global",
            "UDMId:="		, "",
            "MaterialValue:="	, "\"vacuum\"",
            "SurfaceMaterialValue:=", "\"\"",
            "SolveInside:="		, True,
            "ShellElement:="	, False,
            "ShellElementThickness:=", "0mm",
            "ReferenceTemperature:=", "20cel",
            "IsMaterialEditable:="	, True,
            "UseMaterialAppearance:=", False,
            "IsLightweight:="	, False
            ])
        self.obj_names["Rectangle"][0]["Gap"].append(rectangle_name)
        self.obj_names["Rectangle"][1] += 1

        oEditor.Subtract(
            [
                "NAME:Selections",
                "Blank Parts:="		, ",".join([part for part in self.obj_names['Polygon'][0]["Coil_out"]
                                                + self.obj_names['Polygon'][0]["Coil_in"]
                                                if part not in self.subtracted_names]),   #得用“A,B,C"的格式
                "Tool Parts:="		, "Rectangle1"
            ], 
            [
                "NAME:SubtractParameters",
                "KeepOriginals:="	, False,
                "TurnOnNBodyBoolean:="	, True
            ])
        self.subtracted_names.append(rectangle_name)
    
    def create_polylines(self):
        radins1, radouts1, radins2, radouts2 = self.radins1, self.radouts1, self.radins2, self.radouts2
        oEditor = self.ansys_handles[-1]
        x = 4 / np.tan(np.radians(56/58*180 / 2))  #  两边夹角theta = (n-2)/n * 180

        for i in range(len(radins1)-1):
            polyline_name = oEditor.CreatePolyline(
                [
                    "NAME:PolylineParameters",
                    "IsPolylineCovered:="	, True,
                    "IsPolylineClosed:="	, True,
                    [
                        "NAME:PolylinePoints",
                        [
                            "NAME:PLPoint",
                            # "X:="			, f"-{np.sqrt(radins1[i+1]**2-16)}mm",
                            "X:="			, f"-{radins1[i+1]-x}mm",
                            "Y:="			, "-4mm",
                            "Z:="			, "0mm"
                        ],
                        [
                            "NAME:PLPoint",
                            # "X:="			, f"-{np.sqrt(radins1[i]**2-16)}mm",
                            "X:="			, f"-{radins1[i]-x}mm",    
                            "Y:="			, "4mm",
                            "Z:="			, "0mm"
                        ],
                        [
                            "NAME:PLPoint",
                            # "X:="			, f"-{np.sqrt(radouts1[i]**2-16)}mm",   
                            "X:="			, f"-{radouts1[i]-x}mm",  
                            "Y:="			, "4mm",
                            "Z:="			, "0mm"
                        ],
                        [
                            "NAME:PLPoint",
                            # "X:="			, f"-{np.sqrt(radouts1[i+1]**2-16)}mm",
                            "X:="			, f"-{radouts1[i+1]-x}mm",    
                            "Y:="			, "-4mm",
                            "Z:="			, "0mm"
                        ],
                        [
                            "NAME:PLPoint",
                            # "X:="			, f"-{np.sqrt(radins1[i+1]**2-16)}mm",
                            "X:="			, f"-{radins1[i+1]-x}mm",     
                            "Y:="			, "-4mm",
                            "Z:="			, "0mm"
                        ]
                    ],
                    [
                        "NAME:PolylineSegments",
                        [
                            "NAME:PLSegment",
                            "SegmentType:="		, "Line",
                            "StartIndex:="		, 0,
                            "NoOfPoints:="		, 2
                        ],
                        [
                            "NAME:PLSegment",
                            "SegmentType:="		, "Line",
                            "StartIndex:="		, 1,
                            "NoOfPoints:="		, 2
                        ],
                        [
                            "NAME:PLSegment",
                            "SegmentType:="		, "Line",
                            "StartIndex:="		, 2,
                            "NoOfPoints:="		, 2
                        ],
                        [
                            "NAME:PLSegment",
                            "SegmentType:="		, "Line",
                            "StartIndex:="		, 3,
                            "NoOfPoints:="		, 2
                        ]
                    ],
                    [
                        "NAME:PolylineXSection",
                        "XSectionType:="	, "None",
                        "XSectionOrient:="	, "Auto",
                        "XSectionWidth:="	, "0mm",
                        "XSectionTopWidth:="	, "0mm",
                        "XSectionHeight:="	, "0mm",
                        "XSectionNumSegments:="	, "0",
                        "XSectionBendType:="	, "Corner"
                    ]
                ], 
                [
                    "NAME:Attributes",
                    "Name:="		, f"Polyline{self.obj_names['Polyline'][1]}",      # 每生成一个之后会+1，初始是1
                    "Flags:="		, "",
                    "Color:="		, "(143 175 143)",
                    "Transparency:="	, 0,
                    "PartCoordinateSystem:=", "Global",
                    "UDMId:="		, "",
                    "MaterialValue:="	, "\"vacuum\"",
                    "SurfaceMaterialValue:=", "\"\"",
                    "SolveInside:="		, True,
                    "ShellElement:="	, False,
                    "ShellElementThickness:=", "0mm",
                    "ReferenceTemperature:=", "20cel",
                    "IsMaterialEditable:="	, True,
                    "UseMaterialAppearance:=", False,
                    "IsLightweight:="	, False
                ])
            self.obj_names["Polyline"][0]["Connector"].append(polyline_name)
            self.obj_names["Polyline"][1] += 1	# 修改命名       

        for i in range(len(radins2)-1):
            polyline_name = oEditor.CreatePolyline(
                [
                    "NAME:PolylineParameters",
                    "IsPolylineCovered:="	, True,
                    "IsPolylineClosed:="	, True,
                    [
                        "NAME:PolylinePoints",
                        [
                            "NAME:PLPoint",
                            "X:="			, f"-{np.sqrt(radins2[i+1]**2-16)}mm",	    
                            "Y:="			, "-4mm",
                            "Z:="			, "0mm"
                        ],
                        [
                            "NAME:PLPoint",
                            "X:="			, f"-{np.sqrt(radins2[i]**2-16)}mm",    
                            "Y:="			, "4mm",
                            "Z:="			, "0mm"
                        ],
                        [
                            "NAME:PLPoint",
                            "X:="			, f"-{np.sqrt(radouts2[i]**2-16)}mm",     
                            "Y:="			, "4mm",
                            "Z:="			, "0mm"
                        ],
                        [
                            "NAME:PLPoint",
                            "X:="			, f"-{np.sqrt(radouts2[i+1]**2-16)}mm",    
                            "Y:="			, "-4mm",
                            "Z:="			, "0mm"
                        ],
                        [
                            "NAME:PLPoint",
                            "X:="			, f"-{np.sqrt(radins2[i+1]**2-16)}mm",     
                            "Y:="			, "-4mm",
                            "Z:="			, "0mm"
                        ]
                    ],
                    [
                        "NAME:PolylineSegments",
                        [
                            "NAME:PLSegment",
                            "SegmentType:="		, "Line",
                            "StartIndex:="		, 0,
                            "NoOfPoints:="		, 2
                        ],
                        [
                            "NAME:PLSegment",
                            "SegmentType:="		, "Line",
                            "StartIndex:="		, 1,
                            "NoOfPoints:="		, 2
                        ],
                        [
                            "NAME:PLSegment",
                            "SegmentType:="		, "Line",
                            "StartIndex:="		, 2,
                            "NoOfPoints:="		, 2
                        ],
                        [
                            "NAME:PLSegment",
                            "SegmentType:="		, "Line",
                            "StartIndex:="		, 3,
                            "NoOfPoints:="		, 2
                        ]
                    ],
                    [
                        "NAME:PolylineXSection",
                        "XSectionType:="	, "None",
                        "XSectionOrient:="	, "Auto",
                        "XSectionWidth:="	, "0mm",
                        "XSectionTopWidth:="	, "0mm",
                        "XSectionHeight:="	, "0mm",
                        "XSectionNumSegments:="	, "0",
                        "XSectionBendType:="	, "Corner"
                    ]
                ], 
                [
                    "NAME:Attributes",
                    "Name:="		, f"Polyline{self.obj_names['Polyline'][1]}",      # 每生成一个之后会+1，初始是1
                    "Flags:="		, "",
                    "Color:="		, "(143 175 143)",
                    "Transparency:="	, 0,
                    "PartCoordinateSystem:=", "Global",
                    "UDMId:="		, "",
                    "MaterialValue:="	, "\"vacuum\"",
                    "SurfaceMaterialValue:=", "\"\"",
                    "SolveInside:="		, True,
                    "ShellElement:="	, False,
                    "ShellElementThickness:=", "0mm",
                    "ReferenceTemperature:=", "20cel",
                    "IsMaterialEditable:="	, True,
                    "UseMaterialAppearance:=", False,
                    "IsLightweight:="	, False
                ])
            self.obj_names["Polyline"][0]["Connector"].append(polyline_name)
            self.obj_names["Polyline"][1] += 1

    def create_leads_C(self): # 画外线圈导线并且设置边界条件和激励
        oDesign, oEditor = self.ansys_handles[-2], self.ansys_handles[-1]
        in0, out1 = self.radouts2[0], self.radouts1[-1]
        win0, wout1 =self.kwargs["wins"][0], self.kwargs["wouts"][-1]
        x = 4 / np.tan(np.radians(56/58*180 / 2))  #  两边夹角theta = (n-2)/n * 180
        XPositions = [-out1+x, -out1+x, -out1+x+wout1, -in0+x, -in0+x]
        YPositions = [4, 4, 1, -4, -4]
        ZPositions = [0, 2, 2, 2, 0]
        Widths =     [2, wout1, out1-in0-wout1, win0, 2]
        Heights =    [wout1, -5, -2, 5, win0]
        WhichAxiss = ["Y", "Z", "Z", "Z", "Y"]
        for i,(XPosition, YPosition, ZPosition, Width, Height, WhichAxis) in enumerate(zip(XPositions, YPositions, 
                                                                                ZPositions, Widths, 
                                                                                Heights, WhichAxiss)):
            retcangle_name = oEditor.CreateRectangle(
                [
                    "NAME:RectangleParameters",
                    "IsCovered:="		, True,
                    "XStart:="		, f"{XPosition}mm",
                    "YStart:="		, f"{YPosition}mm",
                    "ZStart:="		, f"{ZPosition}mm",
                    "Width:="		, f"{Width}mm",         
                    "Height:="		, f"{Height}mm",        
                    "WhichAxis:="	, WhichAxis
                ], 
                [
                    "NAME:Attributes",
                    "Name:="		, f"Rectangle{self.obj_names['Rectangle'][1]}", # EXITATION
                    "Flags:="		, "",
                    "Color:="		, "(143 175 143)",
                    "Transparency:="	, 0,
                    "PartCoordinateSystem:=", "Global",
                    "UDMId:="		, "",
                    "MaterialValue:="	, "\"vacuum\"",
                    "SurfaceMaterialValue:=", "\"\"",
                    "SolveInside:="		, True,
                    "ShellElement:="	, False,
                    "ShellElementThickness:=", "0mm",
                    "ReferenceTemperature:=", "20cel",
                    "IsMaterialEditable:="	, True,
                    "UseMaterialAppearance:=", False,
                    "IsLightweight:="	, False
                ])
            self.obj_names["Rectangle"][0]["Lead"].append(retcangle_name)   # Rectangle4
            self.obj_names["Rectangle"][1] += 1
        
        oModule = oDesign.GetModule("BoundarySetup")
        # 设置无厚度有限导体边界
        oModule.AssignPerfectE(
            [
                "NAME:PerfE1",
                "Objects:="		, ["Rectangle6","Rectangle5","Rectangle3","Rectangle2"],
                "InfGroundPlane:="	, False
            ])        
        oModule.AssignLumpedRLC(
            [
                "NAME:LumpRLC1",
                "Objects:="		, ["Rectangle4"],
                [
                    "NAME:CurrentLine",
                    "Coordinate System:="	, "Global",
                    "Start:="		, [f"{-out1+x+wout1}mm","0mm","2mm"],
                    "End:="			, [f"{-in0+x}mm","0mm","2mm"]
                ],
                "RLC Type:="		, "Parallel",
                "UseResist:="		, False,
                "UseInduct:="		, False,
                "UseCap:="		, True,
                "Capacitance:="		, "20pF"
            ])
    def create_leads_excitation(self):
        oDesign, oEditor = self.ansys_handles[-2], self.ansys_handles[-1]
        out0, in1 = self.radouts1[0], self.radouts2[-1]
        wout0, win1 = self.kwargs["wouts"][0], self.kwargs["wins"][-1]
        # x = 4 / np.tan(np.radians(56/58*180 / 2))  #  两边夹角theta = (n-2)/n * 180
        XPositions = [-np.sqrt(out0**2-16), -np.sqrt(out0**2-16), -np.sqrt(out0**2-16)+wout0, -np.sqrt(in1**2-16), -np.sqrt(in1**2-16)]
        YPositions = [-4, -4, 1, 4, 4]
        ZPositions = [0, 3, 3, 3, 0]
        Widths =     [3, wout0, np.sqrt(out0**2-16)-np.sqrt(in1**2-16)-wout0, win1, 3]
        Heights =    [wout0, 5, -2, -5, win1]
        WhichAxiss = ["Y", "Z", "Z", "Z", "Y"]
        for i,(XPosition, YPosition, ZPosition, Width, Height, WhichAxis) in enumerate(zip(XPositions, YPositions, 
                                                                                ZPositions, Widths, 
                                                                                Heights, WhichAxiss)):
            retcangle_name = oEditor.CreateRectangle(
                [
                    "NAME:RectangleParameters",
                    "IsCovered:="		, True,
                    "XStart:="		, f"{XPosition}mm",
                    "YStart:="		, f"{YPosition}mm",
                    "ZStart:="		, f"{ZPosition}mm",
                    "Width:="		, f"{Width}mm",         
                    "Height:="		, f"{Height}mm",        
                    "WhichAxis:="	, WhichAxis
                ], 
                [
                    "NAME:Attributes",
                    "Name:="		, f"Rectangle{self.obj_names['Rectangle'][1]}", # EXITATION
                    "Flags:="		, "",
                    "Color:="		, "(143 175 143)",
                    "Transparency:="	, 0,
                    "PartCoordinateSystem:=", "Global",
                    "UDMId:="		, "",
                    "MaterialValue:="	, "\"vacuum\"",
                    "SurfaceMaterialValue:=", "\"\"",
                    "SolveInside:="		, True,
                    "ShellElement:="	, False,
                    "ShellElementThickness:=", "0mm",
                    "ReferenceTemperature:=", "20cel",
                    "IsMaterialEditable:="	, True,
                    "UseMaterialAppearance:=", False,
                    "IsLightweight:="	, False
                ])
            self.obj_names["Rectangle"][0]["Lead"].append(retcangle_name)   # Rectangle9是exitation
            self.obj_names["Rectangle"][1] += 1    

        oModule = oDesign.GetModule("BoundarySetup")
        # 设置无厚度有限导体边界
        oModule.AssignPerfectE(
            [
                "NAME:PerfE2",
                "Objects:="		, ["Rectangle7","Rectangle8","Rectangle10","Rectangle11"],
                "InfGroundPlane:="	, False
            ])     

        # 设置激励
        oModule.AutoIdentifyPorts(
            [
                "NAME:Faces", 
                int(oEditor.GetFaceIDs("Rectangle9")[0])
            ], False, 
            [
                "NAME:ReferenceConductors", 
                "Rectangle10"
            ], "1", True)    

    def create_leads_secondPCB(self):
        oDesign, oEditor = self.ansys_handles[-2], self.ansys_handles[-1]
        out0, in1 = self.radouts1[0], self.radouts2[-1]
        wout0, win1 = self.kwargs["wouts"][0], self.kwargs["wins"][-1]
        in0, out1 = self.radouts2[0], self.radouts1[-1]
        win0, wout1 =self.kwargs["wins"][0], self.kwargs["wouts"][-1]
        x = 4 / np.tan(np.radians(56/58*180 / 2))  #  两边夹角theta = (n-2)/n * 180

        XPositions = [-out1+x, -out1+x, -out1+x+wout1, -in0+x, -in0+x, -np.sqrt(out0**2-16), -np.sqrt(out0**2-16), -np.sqrt(out0**2-16)+wout0, -np.sqrt(in1**2-16), -np.sqrt(in1**2-16)]
        YPositions = [-4, -4, 1, 4, 4, 4, 4, 1, -4, -4]
        ZPositions = [self.PCBspace, self.PCBspace-2, self.PCBspace-2, self.PCBspace-2, self.PCBspace, self.PCBspace, self.PCBspace-3, self.PCBspace-3, self.PCBspace-3, self.PCBspace]
        Widths =     [-2, wout1, out1-in0-wout1, win0, -2, -3, wout0, np.sqrt(out0**2-16)-np.sqrt(in1**2-16)-wout0, win1, -3]
        Heights =    [wout1, 5, -2, -5, win0, wout0, -5, -2, 5, win1]
        WhichAxiss = ["Y", "Z", "Z", "Z", "Y", "Y", "Z", "Z", "Z", "Y"]
        for i,(XPosition, YPosition, ZPosition, Width, Height, WhichAxis) in enumerate(zip(XPositions, YPositions, 
                                                                                ZPositions, Widths, 
                                                                                Heights, WhichAxiss)):
            retcangle_name = oEditor.CreateRectangle(
                [
                    "NAME:RectangleParameters",
                    "IsCovered:="		, True,
                    "XStart:="		, f"{XPosition}mm",
                    "YStart:="		, f"{YPosition}mm",
                    "ZStart:="		, f"{ZPosition}mm",
                    "Width:="		, f"{Width}mm",         
                    "Height:="		, f"{Height}mm",        
                    "WhichAxis:="	, WhichAxis
                ], 
                [
                    "NAME:Attributes",
                    "Name:="		, f"Rectangle{self.obj_names['Rectangle'][1]}", # EXITATION
                    "Flags:="		, "",
                    "Color:="		, "(143 175 143)",
                    "Transparency:="	, 0,
                    "PartCoordinateSystem:=", "Global",
                    "UDMId:="		, "",
                    "MaterialValue:="	, "\"vacuum\"",
                    "SurfaceMaterialValue:=", "\"\"",
                    "SolveInside:="		, True,
                    "ShellElement:="	, False,
                    "ShellElementThickness:=", "0mm",
                    "ReferenceTemperature:=", "20cel",
                    "IsMaterialEditable:="	, True,
                    "UseMaterialAppearance:=", False,
                    "IsLightweight:="	, False
                ])
            self.obj_names["Rectangle"][0]["Lead"].append(retcangle_name)   # Rectangle9是exitation
            self.obj_names["Rectangle"][1] += 1    

        oModule = oDesign.GetModule("BoundarySetup")
        # 设置无厚度有限导体边界
        oModule.AssignPerfectE(
            [
                "NAME:PerfE3",
                "Objects:="		, ["Rectangle12","Rectangle13","Rectangle15","Rectangle16", "Rectangle17", "Rectangle18", "Rectangle20", "Rectangle21"],
                "InfGroundPlane:="	, False
            ])     

        # 设置激励
        oModule.AutoIdentifyPorts(
            [
                "NAME:Faces", 
                int(oEditor.GetFaceIDs("Rectangle19")[0])
            ], False, 
            [
                "NAME:ReferenceConductors", 
                "Rectangle20"
            ], "1", True)    

        oModule.AssignLumpedRLC(
            [
                "NAME:LumpRLC2",
                "Objects:="		, ["Rectangle14"],
                [
                    "NAME:CurrentLine",
                    "Coordinate System:="	, "Global",
                    "Start:="		, [f"{-out1+x+wout1}mm","0mm",f"{self.PCBspace-2}mm"],
                    "End:="			, [f"{-in0+x}mm","0mm",f"{self.PCBspace-2}mm"]
                ],
                "RLC Type:="		, "Parallel",
                "UseResist:="		, False,
                "UseInduct:="		, False,
                "UseCap:="		, True,
                "Capacitance:="		, "20pF"
            ])


    def assign_boundary_singlePCB(self):
        oDesign, oEditor = self.ansys_handles[-2], self.ansys_handles[-1]
        filtered_coil_elements = [elem for elem in self.obj_names["Polygon"][0]["Coil_out"] + self.obj_names["Polygon"][0]["Coil_in"] if elem not in self.subtracted_names] + self.obj_names["Polyline"][0]["Connector"]
        oModule = oDesign.GetModule("BoundarySetup")
        oModule.AssignFiniteCond( 
            [
                "NAME:FiniteCond1",
                "Objects:="		    , filtered_coil_elements,
                "UseMaterial:="		, True,
                "Material:="		, "copper",
                "UseThickness:="	, True,
                "Thickness:="       , "0.09mm",
                # "UseThickness:="	, False,
                "Roughness:="		, "0um",
                "InfGroundPlane:="	, False,
                "IsTwoSided:="		, False,
                "IsInternal:="		, True
            ])

        # 绘制空气盒子
        oEditor.CreateRegion(
            [
                "NAME:RegionParameters",
                "+XPaddingType:="	, "Absolute Position",
                "+XPadding:="		, "500mm",
                "-XPaddingType:="	, "Absolute Position",
                "-XPadding:="		, "-500mm",
                "+YPaddingType:="	, "Absolute Position",
                "+YPadding:="		, "500mm",
                "-YPaddingType:="	, "Absolute Position",
                "-YPadding:="		, "-500mm",
                "+ZPaddingType:="	, "Absolute Position",
                "+ZPadding:="		, "500mm",
                "-ZPaddingType:="	, "Absolute Position",
                "-ZPadding:="		, "-500mm",
                [
                    "NAME:BoxForVirtualObjects",
                    [
                        "NAME:LowPoint", 
                        1, 
                        1, 
                        1
                    ],
                    [
                        "NAME:HighPoint", 
                        -1, 
                        -1, 
                        -1
                    ]
                ]
            ], 
            [
                "NAME:Attributes",
                "Name:="		, "Region",
                "Flags:="		, "Wireframe#",
                "Color:="		, "(143 175 143)",
                "Transparency:="	, 0,
                "PartCoordinateSystem:=", "Global",
                "UDMId:="		, "",
                "MaterialValue:="	, "\"vacuum\"",
                "SurfaceMaterialValue:=", "\"\"",
                "SolveInside:="		, True,
                "ShellElement:="	, False,
                "ShellElementThickness:=", "nan ",
                "ReferenceTemperature:=", "nan ",
                "IsMaterialEditable:="	, True,
                "UseMaterialAppearance:=", False,
                "IsLightweight:="	, False
            ])     
    
    def assign_mesh(self):
        oDesign = self.ansys_handles[-2]
        oModule = oDesign.GetModule("MeshSetup")
        filtered_coil_elements = [elem for elem in self.obj_names["Polygon"][0]["Coil_out"] + self.obj_names["Polygon"][0]["Coil_in"] if elem not in self.subtracted_names]
        all_objects = self.obj_names["Rectangle"][0]["Lead"] + filtered_coil_elements + self.obj_names["Polyline"][0]["Connector"] + self.obj_names["SecondPCB"][0]["Second"] + ["RegularPolyhedron1"]
        # print(all_objects)

        oModule.AssignLengthOp(
            [
                "NAME:Length1",
                "RefineInside:="	, False,
                "Enabled:="		    , True,
                "Objects:="		    , all_objects,
                "RestrictElem:="	, False,
                "NumMaxElem:="		, "1000",
                "RestrictLength:="	, True,
                "MaxLength:="		, "10mm"
            ])

    def analysis_setup_singlepcb(self):
        oDesign = self.ansys_handles[-2]
        oModule = oDesign.GetModule("AnalysisSetup")
        oModule.InsertSetup("HfssDriven", 
            [
                "NAME:Setup1",
                "SolveType:="		, "Broadband",
                [
                    "NAME:MultipleAdaptiveFreqsSetup",
                    "Low:="			, f"{self.frequencyrange[0]}MHz",
                    "High:="		, f"{self.frequencyrange[1]}MHz"
                ],
                "MaxDeltaS:="		, 0.05,
                "MaximumPasses:="	, 5,
                "MinimumPasses:="	, 1,
                "MinimumConvergedPasses:=", 1,
                "PercentRefinement:="	, 30,
                "IsEnabled:="		, True,
                [
                    "NAME:MeshLink",
                    "ImportMesh:="		, False
                ],
                "BasisOrder:="		, 1,
                "DoLambdaRefine:="	, True,
                "DoMaterialLambda:="	, True,
                "SetLambdaTarget:="	, False,
                "Target:="		, 0.3333,
                "UseMaxTetIncrease:="	, False,
                "PortAccuracy:="	, 2,
                "UseABCOnPort:="	, False,
                "SetPortMinMaxTri:="	, False,
                "DrivenSolverType:="	, "Direct Solver",
                "EnhancedLowFreqAccuracy:=", False,
                "SaveRadFieldsOnly:="	, False,
                "SaveAnyFields:="	, True,
                "IESolverType:="	, "Auto",
                "LambdaTargetForIESolver:=", 0.15,
                "UseDefaultLambdaTgtForIESolver:=", True,
                "IE Solver Accuracy:="	, "Balanced",
                "InfiniteSphereSetup:="	, ""
            ])
        oModule.InsertFrequencySweep("Setup1", 
            [
                "NAME:Sweep",
                "IsEnabled:="		, True,
                "RangeType:="		, "LinearCount",
                "RangeStart:="		, f"{self.frequencyrange[0]}MHz",
                "RangeEnd:="		, f"{self.frequencyrange[1]}MHz",
                "RangeCount:="		, 401,
                "Type:="		, "Interpolating",
                "SaveFields:="		, False,
                "SaveRadFields:="	, False,
                "InterpTolerance:="	, 0.5,
                "InterpMaxSolns:="	, 250,
                "InterpMinSolns:="	, 0,
                "InterpMinSubranges:="	, 1,
                "InterpUseS:="		, True,
                "InterpUsePortImped:="	, True,
                "InterpUsePropConst:="	, True,
                "UseDerivativeConvergence:=", False,
                "InterpDerivTolerance:=", 0.2,
                "UseFullBasis:="	, True,
                "EnforcePassivity:="	, True,
                "PassivityErrorTolerance:=", 0.0001,
                "EnforceCausality:="	, False,
                "SMatrixOnlySolveMode:=", "Auto"
            ])
    
    def save_result_singlepcb(self):
        oDesign = self.ansys_handles[-2]
        report = f"index{self.index}-singlePCB"
        report_file =f"{report}.csv"
        self.report_files.append(self.res_path+"/"+report_file)

        oModule = oDesign.GetModule("ReportSetup")
        oModule.CreateReport(report, "Terminal Solution Data", "Rectangular Plot", "Setup1 : Sweep", 
                [
                    "Domain:="		, "Sweep"
                ], 
                [
                    "Freq:="		, ["All"]
                ], 
                [
                    "X Component:="		, "Freq",
                    "Y Component:="		, ["dB20(Zt(Rectangle8_T1,Rectangle8_T1))"]
                ])
        oModule.AddTraces(report, "Setup1 : Sweep", 
            [
                "Domain:="		, "Sweep"
            ], 
            [
                "Freq:="		, ["All"]
            ], 
            [
                "X Component:="		, "Freq",
                "Y Component:="		, ["ang_deg(Zt(Rectangle8_T1,Rectangle8_T1))"]
            ])
        oDesign.Analyze("Setup1 : Sweep")
        oModule = oDesign.GetModule("ReportSetup")
        oModule.ExportToFile(report, self.report_files[-1], False)
    
    def close_design(self):
        oProject = self.ansys_handles[1]
        oProject.Close() 

    def create_secondPCB(self):
        oDesign, oEditor = self.ansys_handles[-2], self.ansys_handles[-1]
        filtered_coil_elements = [elem for elem in self.obj_names["Polygon"][0]["Coil_in"] + self.obj_names["Polygon"][0]["Coil_out"] if elem not in self.subtracted_names]
        oEditor = oDesign.SetActiveEditor("3D Modeler")
        first_pcb = ",".join(self.obj_names["RegularPolyhedron"][0]["Board"]
                            + filtered_coil_elements
                            + self.obj_names["Polyline"][0]["Connector"])
        
        second_pcb = oEditor.DuplicateAroundAxis(
            [
                "NAME:Selections",
                "Selections:="		, first_pcb,
                "NewPartsModelFlag:="	, "Model"
            ], 
            [
                "NAME:DuplicateAroundAxisParameters",
                "CreateNewObjects:="	, True,
                "WhichAxis:="		, "X",
                "AngleStr:="		, "180deg",
                "NumClones:="		, "2"
            ], 
            [
                "NAME:Options",
                "DuplicateAssignments:=", True
            ], 
            [
                "CreateGroupsForNewObjects:=", False
            ])
        self.obj_names["SecondPCB"][0]["Second"] = second_pcb

        oEditor.Move(
            [
                "NAME:Selections",
                "Selections:="		, ",".join(second_pcb),
                "NewPartsModelFlag:="	, "Model"
            ], 
            [
                "NAME:TranslateParameters",
                "TranslateVectorX:="	, "0mm",
                "TranslateVectorY:="	, "0mm",
                "TranslateVectorZ:="	, f"{self.PCBspace}mm"    # EXTENSION
            ])

    def save_results_secondPCB(self):
        oProject, oDesign = self.ansys_handles[1:3]
        report = f"index{self.index}-doublePCB"      
        report_file = f"{report}.csv"    
        self.report_files.append(self.res_path+"/"+report_file)          
        oModule = oDesign.GetModule("ReportSetup")
        oModule.CreateReport(report, "Terminal Solution Data", "Rectangular Plot", "Setup1 : Sweep", 
                [
                    "Domain:="		, "Sweep"
                ], 
                [
                    "Freq:="		, ["All"]
                ], 
                [
                    "X Component:="		, "Freq",
                    "Y Component:="		, ["dB20(Zt(Rectangle8_T1,Rectangle8_T1))"]  
                ])
        oModule.AddTraces(report, "Setup1 : Sweep", 
                [
                    "Domain:="		, "Sweep"
                ], 
                [
                    "Freq:="		, ["All"]
                ], 
                [
                    "X Component:="		, "Freq",
                    "Y Component:="		, ["ang_deg(Zt(Rectangle8_T1,Rectangle8_T1))"]
                ])   
        oProject.Save()
        oDesign.Analyze("Setup1 : Sweep")
        oModule = oDesign.GetModule("ReportSetup")
        oModule.ExportToFile(report, self.report_files[-1], False)

    def run(self):
        self.create_PTFE()
        
        self.create_PCB()
        
        self.create_coils()
        
        self.create_gap()
        
        self.create_polylines()

        self.assign_boundary_singlePCB()
        
        self.create_leads_C()

        self.create_leads_excitation()

        self.analysis_setup_singlepcb()

        self.save_result_singlepcb()
        
        self.create_secondPCB()

        self.create_leads_secondPCB()
        
        self.save_results_secondPCB()
    
    def parse_results(self):
        dB20Z_data = np.loadtxt(self.report_files[2], delimiter=",", skiprows=1)   # 0是粗略阶段，1是精细阶段
        R1_data = np.loadtxt(self.report_files[0], delimiter=",", skiprows=1)
        R2_data = np.loadtxt(self.report_files[1], delimiter=",", skiprows=1) 

        min_indices1 = argrelextrema(dB20Z_data[:, 1], np.less)[0]
        min_indices2 = argrelextrema(dB20Z_data[:, 2], np.less)[0]
        if len(min_indices1) > 0:
            first_valley_index1 = min_indices1[0]
            SRF_out_1 = dB20Z_data[first_valley_index1, 0]
        else:
            idx1 = np.argmin(dB20Z_data[:, 1])
            SRF_out_1 = dB20Z_data[idx1, 0]
        
        if len(min_indices2) > 0:
            first_valley_index2 = min_indices2[0]
            SRF_in_1 = dB20Z_data[first_valley_index2, 0]
        else:
            idx2 = np.argmin(dB20Z_data[:, 2])
            SRF_in_1 = dB20Z_data[idx2, 0] 

        max_indices1 = argrelextrema(dB20Z_data[:, 1], np.greater)[0]
        max_indices2 = argrelextrema(dB20Z_data[:, 2], np.greater)[0]

        if len(max_indices1) > 0:
            first_peek_index1 = max_indices1[0]
            SRF_out_2 = dB20Z_data[first_peek_index1, 0]
        else:
            SRF_out_2 = 11

        if len(max_indices2) > 0:
            first_peek_index2 = max_indices2[0]
            SRF_in_2 = dB20Z_data[first_peek_index2, 0]
        else:
            SRF_in_2 = 11
        
        R1 = R1_data[0, 1]
        R2 = R2_data[0, 1]

        t1 = 5.0
        t2 = 1.0
        t3 = 1.0
        t4 = -1.0
        t5 = -1.0

        obj_value = -t1*(R1+R2)-t2*abs(SRF_out_1-SRF_in_1) - t3*(SRF_out_1+SRF_in_1) - t4*abs(SRF_out_1-SRF_out_2) - t5*abs(SRF_in_1-SRF_in_2)
        self.parsed_results = {"SRFout1": SRF_out_1, "SRFin1": SRF_in_1, "SRFout2":SRF_out_2, "SRFin2":SRF_in_2, "R1": R1, "R2": R2, "obj": obj_value}

        return self.parsed_results
    
def init_ansys(oDesktop, project_name=None):
    
    hfss_app = Hfss(
        projectname=project_name,      # ✅ 指定唯一 Project 名
        designname="HFSSDesign",
        solution_type="DrivenTerminal",
        specified_version="2023.1"
    )

    oDesktop = hfss_app.odesktop
    oDesktop.RestoreWindow()

    oProject = hfss_app.oproject
    oDesign = hfss_app.odesign

    oDesign.SetSolutionType(
        "DrivenTerminal",
        [
            "NAME:Options",
            "EnableAutoOpen:=", False
        ]
    )

    oEditor = oDesign.SetActiveEditor("3D Modeler")

    return oDesktop, oProject, oDesign, oEditor

def run_hfss(config_file, desktop_instance=None, index=None):
    try:
        # 读取配置文件
        with open(config_file, "r+") as f:
            kwargs = json.load(f)
        
        # 初始化 AEDT Desktop 实例
        if desktop_instance is None:
            desktop_instance = Desktop(
                specified_version=f"{kwargs['Ansys_version']}",
                non_graphical=False,
                close_on_exit=False,    # 如果想把ansys关掉的话就设置为True
                student_version=False
            )
        project_name = f"Sample_{index}_HFSS"
        # 初始化 Ansys 环境并解包返回值
        oDesktop, oProject, oDesign, oEditor = init_ansys(desktop_instance.odesktop, project_name=project_name)
        
        # 检查 project_list 是否为空，防止 list index out of range 错误
        kwargs["project_id"] = project_name
        kwargs["project_path"] = desktop_instance.project_path()
        
        # 更新配置文件
        with open(config_file, "w+") as f:
            json.dump(kwargs, f)
        
        # 创建 IPTCoil 对象
        ansys_handles = (oDesktop, oProject, oDesign, oEditor)
        cir_pcb = HPTHFSS(kwargs, ansys_handles, index=index)
        
        # 运行 cir_pcb
        cir_pcb.run()   
        if_success = True
        # cir_pcb.parsed_results = cir_pcb.parse_results() # EXTENSION
    
    except Exception as e:
        # 捕获异常，确保 cir_pcb 存在才赋值错误日志
        print(f"An error occurred: {e}")
        if 'cir_pcb' in locals():
            cir_pcb.error_log = e
        if_success = False
        
    # cir_pcb.close_design()
    # desktop_instance.release_desktop(close_projects=True)

    # 返回 cir_pcb 和 if_success 状态
    return cir_pcb, if_success

