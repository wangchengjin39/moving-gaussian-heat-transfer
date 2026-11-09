# 复现说明

## 1. 运行环境与后处理核对

在仓库根目录打开 PowerShell，再依次运行：

```powershell
python -m pip install -r requirements.txt
```

```powershell
python analysis/compare_probe_speeds.py
```

```powershell
python analysis/validate_cooling_energy.py
```

第一条在已准备好的 Python 环境中安装本项目依赖。第二条读取 `validation` 中两组测点温度，生成速度对比图和汇总表。第三条读取 4 mm/s 组的归档全场温度，核查停止加热后的整板储热。

速度对比看测点，不看整板平均温度。打开 `outputs/probe_1600W_speed_comparison_summary.csv`。可核对的两个峰温是：2 mm/s 为 829.13 ℃，4 mm/s 为 506.78 ℃。时间指标是最高采样时刻。2 mm/s 组有两个并列的最高采样时刻，为 21.853 s 和 21.878 s；4 mm/s 组为 11.103 s。

储热核查看整板。打开 `outputs/cooling_energy_check.json`，比较同一 4 mm/s 组在 20 s 和 40 s 的整板储热，两次都约为 31,997.68 J。报告中的 `change_within_export_rounding_bound` 为 `true`，表示储热差落在导出舍入误差界内。停热后测点温度仍可继续下降，因为热量还在板内重新分布；这不表示整板储热在减少。

复现状态：已在同一台电脑的另一目录复跑。

## 2. Mechanical 中的算例重建

在 Workbench 中建立瞬态热分析。板为长方体：X 为长度 100 mm，Y 为厚度 10 mm，Z 为宽度 40 mm。材料取常物性，比热容 500 J/(kg·K)，热导率 20 W/(m·K)，密度 7800 kg/m³。单元尺寸设为 1.25 mm。热通量加在 Y = 10 mm、面积为 100 mm × 40 mm 的顶面，并且命名这个平面为TOP_HEAT，并发送至求解器；板材初始温度为20℃，边界条件和时间步可以直接引用模型定义。

两组热源都沿同一路径恒速单向移动，求解终止时间都是 40 s。命令按载荷步编号绑定，不要把快组脚本原样复制给慢组。

4 mm/s 组在分析设置中建立两个载荷步，结束时间分别填 20 s 和 40 s。施加顶面热源的命令选择步数 1：

```apdl
/COM,APDL_TEST
ALLSEL,ALL
CSYS,0
CMSEL,S,TOP_HEAT

*GET,ntop,NODE,0,COUNT
*GET,ymin_top,NODE,0,MNLOC,Y
*GET,ymax_top,NODE,0,MXLOC,Y

*VWRITE,ntop,ymin_top,ymax_top
('TOP_CHECK: N=',F10.0,' YMIN=',F12.6,' YMAX=',F12.6)

ALLSEL,ALL
! Moving Gaussian diagnostic: SI units
p_in = 1600
r_heat = 0.005
x_start = 0.010
v_heat = 0.004
zc_heat = 0.020
pi_val = ACOS(-1)
q_peak = p_in/(pi_val*r_heat**2)

! Three independent variables: X, Z, TIME
*DIM,q_move,TABLE,161,65,161,X,Z,TIME

*DO,kk,161,1,-1
    q_move(0,0,kk) = (kk-1)*0.25
*ENDDO

*DO,kk,1,161
    t_heat = (kk-1)*0.25
    xc_heat = x_start+v_heat*t_heat

    *DO,ii,1,161
        q_move(ii,0,kk) = (ii-1)*0.000625
    *ENDDO

    *DO,jj,1,65
        q_move(0,jj,kk) = (jj-1)*0.000625
    *ENDDO

    *DO,ii,1,161
        xx_heat = (ii-1)*0.000625
        *DO,jj,1,65
            zz_heat = (jj-1)*0.000625
            dist_sq = (xx_heat-xc_heat)**2+(zz_heat-zc_heat)**2
            q_move(ii,jj,kk) = q_peak*EXP(-dist_sq/r_heat**2)
        *ENDDO
    *ENDDO
*ENDDO

ALLSEL,ALL
CSYS,0
CMSEL,S,TOP_HEAT
SF,ALL,HFLUX,%q_move%
ALLSEL,ALL
```

撤掉热源的命令选择步数 2，对应 20 s 至 40 s：

```apdl
ALLSEL,ALL
CMSEL,S,TOP_HEAT
SFDELE,ALL,HFLUX
ALLSEL,ALL
```

2 mm/s 组只有一个载荷步，结束时间为 40 s。命令沿用快组第 1 步的脚本，只把 `v_heat` 从 `0.004` 改为 `0.002` m/s：

```apdl
/COM,APDL_TEST
ALLSEL,ALL
CSYS,0
CMSEL,S,TOP_HEAT

*GET,ntop,NODE,0,COUNT
*GET,ymin_top,NODE,0,MNLOC,Y
*GET,ymax_top,NODE,0,MXLOC,Y

*VWRITE,ntop,ymin_top,ymax_top
('TOP_CHECK: N=',F10.0,' YMIN=',F12.6,' YMAX=',F12.6)

ALLSEL,ALL
! Moving Gaussian diagnostic: SI units
p_in = 1600
r_heat = 0.005
x_start = 0.010
v_heat = 0.002
zc_heat = 0.020
pi_val = ACOS(-1)
q_peak = p_in/(pi_val*r_heat**2)

! Three independent variables: X, Z, TIME
*DIM,q_move,TABLE,161,65,161,X,Z,TIME

*DO,kk,161,1,-1
    q_move(0,0,kk) = (kk-1)*0.25
*ENDDO

*DO,kk,1,161
    t_heat = (kk-1)*0.25
    xc_heat = x_start+v_heat*t_heat

    *DO,ii,1,161
        q_move(ii,0,kk) = (ii-1)*0.000625
    *ENDDO

    *DO,jj,1,65
        q_move(0,jj,kk) = (jj-1)*0.000625
    *ENDDO

    *DO,ii,1,161
        xx_heat = (ii-1)*0.000625
        *DO,jj,1,65
            zz_heat = (jj-1)*0.000625
            dist_sq = (xx_heat-xc_heat)**2+(zz_heat-zc_heat)**2
            q_move(ii,jj,kk) = q_peak*EXP(-dist_sq/r_heat**2)
        *ENDDO
    *ENDDO
*ENDDO

ALLSEL,ALL
CSYS,0
CMSEL,S,TOP_HEAT
SF,ALL,HFLUX,%q_move%
ALLSEL,ALL
```

测点设在 (50, 10, 28) mm。求解后导出同一测点的时间–温度文本。后处理结果以 `outputs` 中的对比图和 CSV 为参考，快组再加一项，从求解日志确认第 2 步确实删除了热通量。

## 同机复现结果

快组从头重建，慢组由快组另存后调整。两组均求解至 40 s，求解日志分别记录 0 条错误、3 条警告。快组日志在第 2 个载荷步开始时记录删除了 2560 个热通量面载荷。

| 速度 | 测点采样行数 | 复跑导出与归档对照 |
|---|---:|---|
| 2 mm/s | 1602 | 逐字节一致 |
| 4 mm/s | 1604 | 逐字节一致 |

公开副本中的两份测点导出与归档逐字节一致，其 SHA-256 见[复现数据清单](../validation/source_manifest.json)。求解器输入和模型元数据仅移除了本机路径。完整求解日志保存在本地归档。
