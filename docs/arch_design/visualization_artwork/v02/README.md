# OmicsWeave 绘图艺术素材 v0.2

本轮是图标、插画和界面样式验证。三套图标沿用同一批语义 ID，28 个中有 24 个是新增语义，4 个是原有概念的重绘稿。素材独立于当前 Python 绘图主题和摘要渲染器，便于比较版本和独立引用。

## 设计取舍

| 场景 | 本次采用的设计 |
| --- | --- |
| 论文机制图 | 纯白版使用浅色填充和深色轮廓。按图标类型配色，保留独立图例、效应色标和证据等级符号。网格用于布局或低密度示意图。 |
| 科学手绘 | 在暖灰、陶土、矿物绿和淡紫上绘制有机轮廓、复描线和局部排线。技术演示板采用深底、大标题、单一重点引线及放大的细胞器局部。 |
| 像素概念与界面 | 用实验室地面、格栅、阴影、模块化管道和嵌入式样本槽构建场景层级。32×32 图标用整数网格和硬边轮廓表达不同材料。像素字用于标题、图层和状态标签。 |

图标色板包含水青、陶土、琥珀、淡紫、天蓝和灰绿六个色族，每个色族定义四级明度。分支、可逆、跨区室交换和证据冲突主要以形状区分。量值和统计含义应由研究结果及其图例指定。代谢物符号为概念图标，具体化学结构应由独立结构图表示。

## 交付文件

- [资产清单](manifest.json) 记录图标 ID、语义类别、所属色族、视口、各风格路径和 UI 管件接点。
- [三风格图标目录](previews/icon_catalog.svg) 对照展示 84 个 SVG。像素版同时提供 28 个 32×32 透明 PNG。
- [像素 UI 样板](previews/pixel_ui.svg) 展示组学图层、示踪样本和区室示意。8 个组件包含实验室地板、格栅、直管、弯头、三通、仪表框、图层控件和样本槽，每个均附 SVG 与原尺寸 PNG。管件的 W/E/N/S 接点按组件画布描述。样板以静态预览呈现。
- [科学手绘演示板](previews/ink_keynote.svg) 使用两个可独立放置的透明底 SVG，[肝细胞示意](components/ink/hepatocyte_cutaway.svg) 和 [线粒体局部](components/ink/mitochondrial_inset.svg)。
- 像素 UI 使用方舟像素 12 px 比例版的拉丁与简体中文字形。预览中的排版文字保存为矢量轮廓，可在离线环境打开。随素材提供上游原始拉丁 WOFF2 和 [OFL-1.1 许可文本](font/OFL.txt)。简体中文字体可从同一官方发行包提取。 [TakWolf 官方发行页](https://github.com/TakWolf/ark-pixel-font/releases/tag/2026.11.10)

## 参考来源与使用方式

- Smallville 和 LimeZu 参考地面材质、道具分层与模块化场景组织。图标、地砖和界面在本项目内重新绘制，源素材包和项目截图未随交付分发。 [Stanford HAI](https://hai.stanford.edu/news/computational-agents-exhibit-believable-humanlike-behavior) · [LimeZu Modern Interiors](https://limezu.itch.io/moderninteriors)
- Factorio 参考直管、转角和分支的可拼接组织方式。本次原创管件用于示意界面中的连线，实际反应、流向和统计证据由 OmicsWeave 数据及图例定义。 [Factorio 管道说明](https://wiki.factorio.com/Pipe)
- NVIDIA GTC 参考深色舞台、大标题、亮色引导线与局部放大的视觉层次。手绘轮廓、配色组合和版式均在本项目内制作。 [NVIDIA GTC Keynote](https://blogs.nvidia.com/blog/nvidia-keynote-at-gtc-2025-ai-news-live-updates/)

两个组合样板均为无实验数据的设计示意。艺术素材采用 AI 辅助编码。选用资产进入稿件前，应结合正式样图的刊物要求、终稿尺寸和可读性复核。项目当前的 SVG 加载器只支持现有资源子集，后续接入应明确 SVG 导入或转换路径。
