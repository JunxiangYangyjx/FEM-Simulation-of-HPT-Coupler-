import numpy as np
import json
from scipy.signal import argrelextrema
from pyaedt import Maxwell3d, Desktop

'''
1. 先建模外线圈五匝，然后建模内线圈五匝（平面）
2. 创建gap，用线组成平面连接内外匝
3. 创建新材料PTFE
4. 然后建模PCB基板（平面）
5. 先把平面线圈变成体
6. 复制平面线圈
7. 把PCB基板变成体
8. 复制PCB基板
9. 创建导线使得线圈闭合
10. 设置材料给PCB基板和线圈
11. 切分平面，删除掉其中一个平面，加入电流激励，重复四次
12. 创建空气域（用region还是添加边界条件）
13. 创建矩阵，四个电流激励
14. InsertSetup
15. 创建报告
16. 导出结果
'''

class HPT(object):
    def __init__(self, kwargs, ansys_handles, index=None):
        super(HPT, self).__init__()
        self.kwargs = kwargs	# 参数
        self.index = index		# 第几次仿真
        self.obj_names = {"RegularPolyhedron":[{"Board": []}, 1],
                          "Polygon":[{"Coil_out":[], "Coil_in":[]},1],
                          "Rectangle":[{"Gap":[]},1],
                          "Box":[{"Lead": []},1],
                          "Polyline": [{"Connector": []}, 1],
                        #   "Exitation":[{"Exitation":[]},2],
                          "SecondPCB":[{"Second":[]},1]
                          } # EXTENSION
        self.subtracted_names = []
        self.boundary = []
        self.PCBspace = 10 # mm
        self.roughfrequencyrange = [1, 11]  # 粗略扫频的频率范围, 要看看SRF在什么地方，如果比较小的话感觉不太需要扫描两次
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

    def create_leads_out(self): # 画外线圈导线并且设置边界条件和激励
        oEditor = self.ansys_handles[-1]
        out0, out1 = self.radouts1[0], self.radouts1[-1]
        wout0, wout1 =self.kwargs["wouts"][0], self.kwargs["wouts"][-1]
        x = 4 / np.tan(np.radians(56/58*180 / 2))  #  两边夹角theta = (n-2)/n * 180
        XPositions = [-out0+x, -out0+x, -out0+x+wout0, -out1+x, -out1+x, -out0+x, -out0+x, -out0+x+wout0, -out1+x, -out1+x]
        YPositions = [-4, -4, 1, 4, 4, 4, 4, 1, -4, -4]
        ZPositions = [0, 2, 2, 2, 0, self.PCBspace, self.PCBspace-2, self.PCBspace-2, self.PCBspace-2, self.PCBspace]
        XSize =      [wout0, wout0, out0-out1-wout0, wout1, wout1, wout0, wout0, out0-out1-wout0, wout1, wout1]
        YSize =      [0.1, 5, -2, -5, -0.1, -0.1, -5, -2, 5, 0.1]
        ZSize =      [2, 0.1, 0.1, 0.1, 2, -2, -0.1, -0.1, -0.1, -2]
        for i,(XPosition, YPosition, ZPosition, XSize, YSize, ZSize) in enumerate(zip(XPositions, YPositions, 
                                                                                ZPositions, XSize, 
                                                                                YSize, ZSize)):
            box_name = oEditor.CreateBox(
                [
                    "NAME:BoxParameters",
                    "XPosition:="		, f"{XPosition}mm",
                    "YPosition:="		, f"{YPosition}mm",
                    "ZPosition:="		, f"{ZPosition}mm",
                    "XSize:="		    , f"{XSize}mm",         
                    "YSize:="		    , f"{YSize}mm",        
                    "ZSize:="	        , f"{ZSize}mm",  
                ], 
                [
                    "NAME:Attributes",
                    "Name:="		, f"Box{self.obj_names['Box'][1]}", # EXITATION
                    "Flags:="		, "",
                    "Color:="		, "(143 175 143)",
                    "Transparency:="	, 0,
                    "PartCoordinateSystem:=", "Global",
                    "UDMId:="		, "",
                    "MaterialValue:="	, "\"copper\"",
                    "SurfaceMaterialValue:=", "\"\"",
                    "SolveInside:="		, True,
                    "ShellElement:="	, False,
                    "ShellElementThickness:=", "0mm",
                    "ReferenceTemperature:=", "20cel",
                    "IsMaterialEditable:="	, True,
                    "UseMaterialAppearance:=", False,
                    "IsLightweight:="	, False
                ])
            self.obj_names["Box"][0]["Lead"].append(box_name)  
            self.obj_names["Box"][1] += 1

    def create_leads_in(self): # 删除外线圈激励，画内线圈导线并设置导线和激励
        oEditor = self.ansys_handles[-1]
        in0, in1 = self.radouts2[0], self.radouts2[-1]
        win0, win1 = self.kwargs["wins"][0], self.kwargs["wins"][-1]
        # x = 4 / np.tan(np.radians(56/58*180 / 2))  #  两边夹角theta = (n-2)/n * 180
        XPositions = [-np.sqrt(in0**2-16), -np.sqrt(in0**2-16), -np.sqrt(in0**2-16)+win0, -np.sqrt(in1**2-16), -np.sqrt(in1**2-16), -np.sqrt(in0**2-16), -np.sqrt(in0**2-16), -np.sqrt(in0**2-16)+win0, -np.sqrt(in1**2-16), -np.sqrt(in1**2-16)]
        YPositions = [-4, -4, 1, 4, 4, 4, 4, 1, -4, -4]
        ZPositions = [0, 2, 2, 2, 0, self.PCBspace, self.PCBspace-2, self.PCBspace-2, self.PCBspace-2, self.PCBspace]
        XSize =      [win0, win0, np.sqrt(in0**2-16)-np.sqrt(in1**2-16)-win0, win1, win1, win0, win0, np.sqrt(in0**2-16)-np.sqrt(in1**2-16)-win0, win1, win1]
        YSize =      [0.1, 5, -2, -5, -0.1, -0.1, -5, -2, 5, 0.1]
        ZSize =      [2, 0.1, 0.1, 0.1, 2, -2, -0.1, -0.1, -0.1, -2]
        for i,(XPosition, YPosition, ZPosition, XSize, YSize, ZSize) in enumerate(zip(XPositions, YPositions, 
                                                                                ZPositions, XSize, 
                                                                                YSize, ZSize)):
            box_name = oEditor.CreateBox(
                [
                    "NAME:BoxParameters",
                    "XPosition:="		, f"{XPosition}mm",
                    "YPosition:="		, f"{YPosition}mm",
                    "ZPosition:="		, f"{ZPosition}mm",
                    "XSize:="		    , f"{XSize}mm",         
                    "YSize:="		    , f"{YSize}mm",        
                    "ZSize:="	        , f"{ZSize}mm",  
                ], 
                [
                    "NAME:Attributes",
                    "Name:="		, f"Box{self.obj_names['Box'][1]}", # EXITATION
                    "Flags:="		, "",
                    "Color:="		, "(143 175 143)",
                    "Transparency:="	, 0,
                    "PartCoordinateSystem:=", "Global",
                    "UDMId:="		, "",
                    "MaterialValue:="	, "\"copper\"",
                    "SurfaceMaterialValue:=", "\"\"",
                    "SolveInside:="		, True,
                    "ShellElement:="	, False,
                    "ShellElementThickness:=", "0mm",
                    "ReferenceTemperature:=", "20cel",
                    "IsMaterialEditable:="	, True,
                    "UseMaterialAppearance:=", False,
                    "IsLightweight:="	, False
                ])
            self.obj_names["Box"][0]["Lead"].append(box_name)  
            self.obj_names["Box"][1] += 1 

    def create_region(self):
        oEditor = self.ansys_handles[-1]
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

        
    def sweep_copper(self):
        oEditor = self.ansys_handles[-1]
        filtered_coil_elements = [elem for elem in self.obj_names["Polygon"][0]["Coil_in"] + self.obj_names["Polygon"][0]["Coil_out"] if elem not in self.subtracted_names]
        first_pcb_coil = ",".join(filtered_coil_elements + self.obj_names["Polyline"][0]["Connector"])        
        second_pcb_coil = ",".join([name + "_1" for name in first_pcb_coil.split(",")])
        oEditor.SweepAlongVector(
            [
                "NAME:Selections",
                "Selections:="		, first_pcb_coil,
                "NewPartsModelFlag:="	, "Model"
            ], 
            [
                "NAME:VectorSweepParameters",
                "DraftAngle:="		, "0deg",
                "DraftType:="		, "Round",
                "CheckFaceFaceIntersection:=", False,
                "ClearAllIDs:="		, False,
                "SweepVectorX:="	, "0mm",
                "SweepVectorY:="	, "0mm",
                "SweepVectorZ:="	, "0.035mm"
            ])
        oEditor.SweepAlongVector(
            [
                "NAME:Selections",
                "Selections:="		, second_pcb_coil,
                "NewPartsModelFlag:="	, "Model"
            ], 
            [
                "NAME:VectorSweepParameters",
                "DraftAngle:="		, "0deg",
                "DraftType:="		, "Round",
                "CheckFaceFaceIntersection:=", False,
                "ClearAllIDs:="		, False,
                "SweepVectorX:="	, "0mm",
                "SweepVectorY:="	, "0mm",
                "SweepVectorZ:="	, "-0.035mm"
            ])       
        oEditor.AssignMaterial(
            [
                "NAME:Selections",
                "AllowRegionDependentPartSelectionForPMLCreation:=", True,
                "AllowRegionSelectionForPMLCreation:=", True,
                "Selections:="		, first_pcb_coil,
            ], 
            [
                "NAME:Attributes",
                "MaterialValue:="	, "\"copper\"",
                "SolveInside:="		, True,
                "ShellElement:="	, False,
                "ShellElementThickness:=", "nan ",
                "ReferenceTemperature:=", "nan ",
                "IsMaterialEditable:="	, True,
                "UseMaterialAppearance:=", False,
                "IsLightweight:="	, False
            ])
        oEditor.AssignMaterial(
            [
                "NAME:Selections",
                "AllowRegionDependentPartSelectionForPMLCreation:=", True,
                "AllowRegionSelectionForPMLCreation:=", True,
                "Selections:="		, second_pcb_coil,
            ], 
            [
                "NAME:Attributes",
                "MaterialValue:="	, "\"copper\"",
                "SolveInside:="		, True,
                "ShellElement:="	, False,
                "ShellElementThickness:=", "nan ",
                "ReferenceTemperature:=", "nan ",
                "IsMaterialEditable:="	, True,
                "UseMaterialAppearance:=", False,
                "IsLightweight:="	, False
            ])

    def assign_excitation(self):
        oDesign, oEditor = self.ansys_handles[-2], self.ansys_handles[-1]
        inner_coil_firstpcb = self.obj_names["Polygon"][0]["Coil_in"][0]
        inner_coil_secondpcb = inner_coil_firstpcb + "_1"
        all_inner_coil = [inner_coil_firstpcb, inner_coil_secondpcb]
        inner_coil_section = [
            f"{coil}_Section1"
            for coil in all_inner_coil
        ]
        inner_coil_section_separate = [
            f"{coil}_Section1_Separate1"
            for coil in all_inner_coil
        ]
        oEditor.Section(
            [
                "NAME:Selections",
                "Selections:="		, "Polygon1_1,Polygon1",
                "NewPartsModelFlag:="	, "Model"
            ], 
            [
                "NAME:SectionToParameters",
                "CreateNewObjects:="	, True,
                "SectionPlane:="	, "YZ",
                "SectionCrossObject:="	, False
            ])
        oEditor.Section(
            [
                "NAME:Selections",
                "Selections:="		, ",".join(all_inner_coil),
                "NewPartsModelFlag:="	, "Model"
            ], 
            [
                "NAME:SectionToParameters",
                "CreateNewObjects:="	, True,
                "SectionPlane:="	, "YZ",
                "SectionCrossObject:="	, False
            ])
        
        oEditor.SeparateBody(
            [
                "NAME:Selections",
                "Selections:="		, "Polygon1_1_Section1,Polygon1_Section1",
                "NewPartsModelFlag:="	, "Model"
            ], 
            [
                "CreateGroupsForNewObjects:=", False
            ])
        oEditor.SeparateBody(
            [
                "NAME:Selections",
                "Selections:="		, ",".join(inner_coil_section),
                "NewPartsModelFlag:="	, "Model"
            ], 
            [
                "CreateGroupsForNewObjects:=", False
            ])        
        oEditor.Delete(
            [
                "NAME:Selections",
                "Selections:="		, "Polygon1_1_Section1_Separate1,Polygon1_Section1_Separate1"
            ])
        oEditor.Delete(
            [
                "NAME:Selections",
                "Selections:="		, ",".join(inner_coil_section_separate),
            ])
        oModule = oDesign.GetModule("BoundarySetup")        
        oModule.AssignCurrent(
            [
                "NAME:Current1WAI",
                "Objects:="		, ["Polygon1_Section1"],
                "Phase:="		, "0deg",
                "Current:="		, "1A",
                "IsSolid:="		, True,
                "Point out of terminal:=", False
            ])       
        oModule.AssignCurrent(
            [
                "NAME:Current2WAI",
                "Objects:="		, ["Polygon1_1_Section1"],
                "Phase:="		, "0deg",
                "Current:="		, "1A",
                "IsSolid:="		, True,
                "Point out of terminal:=", True
            ])
        oModule.AssignCurrent(
            [
                "NAME:Current1NEI",
                "Objects:="		, [inner_coil_section[0]],
                "Phase:="		, "0deg",
                "Current:="		, "1A",
                "IsSolid:="		, True,
                "Point out of terminal:=", False
            ])
        oModule.AssignCurrent(
            [
                "NAME:Current2NEI",
                "Objects:="		, [inner_coil_section[1]],
                "Phase:="		, "0deg",
                "Current:="		, "1A",
                "IsSolid:="		, True,
                "Point out of terminal:=", True
            ])
        oModule = oDesign.GetModule("MaxwellParameterSetup")
        oModule.AssignMatrix(
            [
                "NAME:Matrix1",
                [
                    "NAME:MatrixEntry",
                    [
                        "NAME:MatrixEntry",
                        "Source:="		, "Current1NEI",
                        "NumberOfTurns:="	, "1"
                    ],
                    [
                        "NAME:MatrixEntry",
                        "Source:="		, "Current1WAI",
                        "NumberOfTurns:="	, "1"
                    ],
                    [
                        "NAME:MatrixEntry",
                        "Source:="		, "Current2NEI",
                        "NumberOfTurns:="	, "1"
                    ],
                    [
                        "NAME:MatrixEntry",
                        "Source:="		, "Current2WAI",
                        "NumberOfTurns:="	, "1"
                    ]
                ],
                [
                    "NAME:MatrixGroup"
                ]
            ])
        
    def analysis_setup(self):
        oDesign = self.ansys_handles[-2]
        oModule = oDesign.GetModule("AnalysisSetup")
        oModule.InsertSetup("EddyCurrent", 
            [
                "NAME:Setup1",
                "Enabled:="		, True,
                [
                    "NAME:MeshLink",
                    "ImportMesh:="		, False
                ],
                "MaximumPasses:="	, 5,
                "MinimumPasses:="	, 1,
                "MinimumConvergedPasses:=", 1,
                "PercentRefinement:="	, 30,
                "SolveFieldOnly:="	, False,
                "PercentError:="	, 1,
                "SolveMatrixAtLast:="	, True,
                "UseNonLinearIterNum:="	, False,
                "CacheSaveKind:="	, "Delta",
                "ConstantDelta:="	, "0s",
                "UseIterativeSolver:="	, False,
                "RelativeResidual:="	, 1E-05,
                "NonLinearResidual:="	, 0.0001,
                "SmoothBHCurve:="	, False,
                "Frequency:="		, "20Hz",
                "HasSweepSetup:="	, False,
                "UseHighOrderShapeFunc:=", False,
                "UseMuLink:="		, False
            ])

    def save_results(self):
        oProject, oDesign = self.ansys_handles[1:3]
        oModule = oDesign.GetModule("ReportSetup")

        # ===== 1. 生成唯一 report 名和文件名 =====
        report = f"index{self.index}_L"
        report_file = f"{report}.csv"
        save_path = self.res_path + "/" + report_file
        self.report_files.append(save_path)

        # ===== 2. 创建 L Matrix Report =====
        oModule.CreateReport(
            report, "EddyCurrent", "Data Table", "Setup1 : LastAdaptive", [],
            ["Freq:=", ["All"]],
            ["X Component:=", "Freq",
            "Y Component:=", ["Matrix1.L(Current1NEI,Current1NEI)"]]
        )

        oModule.AddTraces(report, "Setup1 : LastAdaptive", [],
            ["Freq:=", ["All"]],
            ["X Component:=", "Freq",
            "Y Component:=", ["Matrix1.L(Current1WAI,Current1WAI)"]]
        )

        oModule.AddTraces(report, "Setup1 : LastAdaptive", [],
            ["Freq:=", ["All"]],
            ["X Component:=", "Freq",
            "Y Component:=", ["Matrix1.L(Current1WAI,Current1NEI)"]]
        )

        oModule.AddTraces(report, "Setup1 : LastAdaptive", [],
            ["Freq:=", ["All"]],
            ["X Component:=", "Freq",
            "Y Component:=", ["Matrix1.L(Current1NEI,Current2NEI)"]]
        )

        oModule.AddTraces(report, "Setup1 : LastAdaptive", [],
            ["Freq:=", ["All"]],
            ["X Component:=", "Freq",
            "Y Component:=", ["Matrix1.L(Current1NEI,Current2WAI)"]]
        )

        oModule.AddTraces(report, "Setup1 : LastAdaptive", [],
            ["Freq:=", ["All"]],
            ["X Component:=", "Freq",
            "Y Component:=", ["Matrix1.L(Current1WAI,Current2WAI)"]]
        )

        # ===== 3. Analyze + Export =====
        oProject.Save()
        oDesign.Analyze("Setup1")

        oModule = oDesign.GetModule("ReportSetup")
        oModule.ExportToFile(report, save_path, False)


    def run(self):
        self.create_PTFE()
        
        self.create_PCB()
        
        self.create_coils()
        
        self.create_gap()
        
        self.create_polylines()

        self.create_secondPCB()
        
        self.sweep_copper()
        
        self.create_leads_out()
        
        self.create_leads_in()

        self.create_region()
        
        self.assign_excitation()
        
        self.analysis_setup()
        
        self.save_results()

    
def init_ansys(oDesktop,
               solution_type="EddyCurrent",
               project_name=None,
               design_name="Design1"):

    maxwell_app = Maxwell3d(
        projectname=project_name,   # ✅ 关键修改
        designname=design_name,
        solution_type=solution_type,
        specified_version="2023.1"
    )

    oDesktop = maxwell_app.odesktop
    oDesktop.RestoreWindow()

    oProject = maxwell_app.oproject
    oDesign = maxwell_app.odesign

    oEditor = oDesign.SetActiveEditor("3D Modeler")

    return oDesktop, oProject, oDesign, oEditor

def run_maxwell(config_file, desktop_instance=None, index=None):
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
        
        # 初始化 Ansys 环境并解包返回值
        project_name = f"Sample_{index}_Maxwell"

        oDesktop, oProject, oDesign, oEditor = init_ansys(
            desktop_instance.odesktop,
            solution_type="EddyCurrent",
            project_name=project_name,
            design_name="MaxwellDesign"
        )       
        
        # 检查 project_list 是否为空，防止 list index out of range 错误
        project_list = desktop_instance.project_list()
        if project_list:
            kwargs["project_id"] = project_list[-1]
        else:
            raise ValueError("No project found in the desktop instance.")
        
        # 记录项目路径
        kwargs["project_path"] = desktop_instance.project_path()
        
        # 更新配置文件
        with open(config_file, "w+") as f:
            json.dump(kwargs, f)
        
        # 创建 IPTCoil 对象
        ansys_handles = (oDesktop, oProject, oDesign, oEditor)
        cir_pcb = HPT(kwargs, ansys_handles, index=index)
        
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

