$WORKDIR="M:\MoAve\Dublin\Area2"
$IMGDIR="$WORKDIR\images_undist"
$RESDIR="$WORKDIR\res"
$SIGMA=1

$INPOSE="$IMGDIR\pose.qinv2"
$OUTPOSE="$IMGDIR\pose_mvs_sol.qinv2"

$GRAPH_NAME="graph_reproject.txt"
$GRAPH_SOV_NAME="graph_reproject_sol.txt"


$BINDIR="C:\ssdev\method\motionaverage\x64\Release"
$SCRIPTDIR="C:\ssdev\method\motionaverage\scripts"

;; ## Create translate graph Build graph and estimate offset from MSP project
& $BINDIR\create_translate_graph_MSP.exe $RESDIR $INPOSE $GRAPH_NAME -v 0 --minrays 3 --sigma $SIGMA

;; ## Minimize graph
& $BINDIR\solve_translate_graph.exe $RESDIR\$GRAPH_NAME $RESDIR\$GRAPH_SOV_NAME

& $BINDIR\apply_translate_graph_MSP.exe $RESDIR $GRAPH_SOV_NAME --inqinv2 $INPOSE --outqinv2 $OUTPOSE


;; ############### Densify
$MSP2COLMAP="J:\Temp_for_SS\msp2colmap\bin\msp2colmap.exe"
;; ########### Convert to COLMAP
$COLMAP_FOLDER="$RESDIR\colmap"
$COLMAP_BEFORE_FOLDER="$COLMAP_FOLDER\before"
$COLMAP_AFTER_FOLDER="$COLMAP_FOLDER\after"
& $MSP2COLMAP -i $INPOSE -a $OUTPOSE -o $RESDIR -s 1000

;; ########### Convert to COLMAP
$COLMAP2MSP="C:\ODM\SuperBuild\install\bin\InterfaceCOLMAP.exe"

mkdir "$COLMAP_BEFORE_FOLDER\sparse"
xcopy "$COLMAP_BEFORE_FOLDER\*.txt" "$COLMAP_BEFORE_FOLDER\sparse"
mkdir "$COLMAP_AFTER_FOLDER\sparse"
copy "$COLMAP_AFTER_FOLDER\*.txt" "$COLMAP_AFTER_FOLDER\sparse"

& $COLMAP2MSP -w "$COLMAP_BEFORE_FOLDER" --image-folder "$IMGDIR" -i "$COLMAP_BEFORE_FOLDER" -o "scene.mvs"
& $COLMAP2MSP -w "$COLMAP_AFTER_FOLDER" --image-folder "$IMGDIR" -i "$COLMAP_AFTER_FOLDER" -o "scene.mvs"
$MVSDENSIFY="C:\ODM\SuperBuild\install\bin\DensifyPointCloud.exe"
& $MVSDENSIFY -w "$COLMAP_BEFORE_FOLDER" -i scene.mvs --resolution-level 0 -v 0
& $MVSDENSIFY -w "$COLMAP_AFTER_FOLDER" -i scene.mvs --resolution-level 0 -v 0

#### Registration
$CLOUD_BEFIRE="$COLMAP_BEFORE_FOLDER\scene_dense.ply"
$CLOUD_AFTER="$COLMAP_AFTER_FOLDER\scene_dense.ply"
$GT_LIDAR="$WORKDIR\lidar_ROI.las"

