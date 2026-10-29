"""MFA 求解结果适配器：旧版自研软件导出格式 → 标准格式。"""

from .legacy import LegacyConversion, convert_legacy_mfa, legacy_flux_matrix, legacy_mid_matrices

__all__ = ["LegacyConversion", "convert_legacy_mfa", "legacy_flux_matrix", "legacy_mid_matrices"]
