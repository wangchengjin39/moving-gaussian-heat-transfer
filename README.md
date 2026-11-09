# 移动高斯热源下的板材导热

本项目比较同一块板、同一个测点在两种焊接速度下的温度响应。两组都使用 1600 W 移动高斯热源，沿同一条路径移动；2 mm/s 加热 40 s，4 mm/s 加热 20 s，观测都到 40 s。主要发现是：速度更低时测点峰值更高、出现更晚，2 mm/s 的峰温为 829.13 ℃，4 mm/s 为 506.78 ℃。

原始数据来自 Ansys Mechanical 2026 R1 Student，后处理使用 Python。材料是常物性练习材料，不是某个已验证牌号。本项目不是焊接质量预测模型。

## 复现状态

Python 后处理已在同一台电脑的另一个目录复跑。在仓库根目录运行下面两条脚本，即可由 `validation` 中的归档导出重新生成对比图、汇总表和储热核查。

Mechanical 的复现在同一台电脑完成，分为两组。4 mm/s 组已从头重建。2 mm/s 组由该快组另存后调整工况。两组测点导出分别与归档逐字节一致。

## 快速运行

已在 Windows、Python 3.14.4 及 `requirements.txt` 列出的版本中测试。先创建自己的 Python 环境，再在仓库根目录执行：

```text
python -m pip install -r requirements.txt
python analysis/compare_probe_speeds.py
python analysis/validate_cooling_energy.py
```

脚本从项目目录读取随附数据，无需安装 Ansys。输出写入 `outputs`；`validation` 中的原始证据不被重写。若缺少依赖，先完成安装。目前未在全新电脑、全新环境中验证安装。

## 主要结果

| 速度/mm·s⁻¹ | 加热终止/s | 测点峰温/℃ | 最高采样时刻/s | 40s测点温度/℃ |
|---|---:|---:|---|---:|
| 2 | 40 | 829.13 | 21.853、21.878（舍入后并列） | 539.02 |
| 4 | 20 | 506.78 | 11.103 | 271.56 |

两组均观察至 40 s，测点为全局 (50, 10, 28) mm，距焊道中心线 8 mm。峰值采样时刻不能当作连续精确峰时或严格峰时区间。4 mm/s 组 20 s 与 40 s 的体积平均温度均约 225.113 ℃，储热均约 31997.68 J；差在导出舍入误差界内。

## 文件导航

- [模型定义](docs/model_definition.md)：几何、材料、热源、测点、单位和求解设置。
- [复现步骤](docs/reproduction.md)：Python 后处理命令，以及 Mechanical 中重建两组工况的步骤。
- [技术报告](docs/technical_report.md)：研究问题、结果、储热核查、局限和个人工作。
- `analysis/compare_probe_speeds.py`：读取归档探针数据，核对共同网格、材料和测点，输出对比图、CSV 和 `outputs/speed_comparison_check.json`。
- `analysis/validate_cooling_energy.py`：按实际二十节点单元连接积分两个时刻的储热，输出 `outputs/cooling_energy_check.json`。
- [复现数据与哈希清单](validation/source_manifest.json)包含脱敏的求解器输入与模型元数据、测点和节点温度导出，以及全场温度极值历程。
- `docs/speed_comparison.md`：速度分析的订正稿。
- [储热核查](docs/cooling_energy_check.md)说明积分方法与结果。

## 模型局限

模型为 100 mm × 40 mm × 10 mm 板，Y 为厚度，初温 20 ℃。材料取常数 k = 20 W/(m·K)、ρ = 7800 kg/m³、c = 500 J/(kg·K)，外表面绝热，无相变。高温区的常物性和无相变假设限制定量解释，不能把约 3143 ℃ 的全场数值当作真实焊接温度。目前没有实验验证。

4 mm/s 组在 20 s 删除热通量，核对结果见[复现步骤](docs/reproduction.md)。公开版的求解器输入和模型元数据已移除本机路径，供 Python 脚本核对模型与网格。在 Mechanical 中重建工况请按复现步骤操作。

低功率网格、时间步和 TABLE 采样已有有限比较，完整证据仍在个人项目归档，这里不声称全场收敛或焊接质量提升。仿射二十节点形函数体积分与节点温度算术平均不同；不能把极小储热差的百分比包装成求解器精度。

## 工作归属

建模、网格和测点设置、参数求解与导出、理论预判和核心解释由项目学习者完成，施热和撤热的载荷步曾在指导下修正。重复比较、积分脚本和 README 初稿由助手整理，项目学习者完成了交付版修订。移动热源的 APDL 命令来自助手提供的脚本，学习者按组别修改了速度等参数。
