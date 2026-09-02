# Binarios de AutoDock Vina — procedencia y checksums

Los binarios no se versionan en git (`tools/vina/*.exe` está en `.gitignore`); este archivo sí, para que la descarga sea reproducible.

**Fuente:** release oficial `ccsb-scripps/AutoDock-Vina` v1.2.7 (GitHub).
**Fecha de descarga:** 2026-08-21.

| Archivo | Tamaño (bytes) | SHA-256 | URL |
|---|---|---|---|
| `vina_1.2.7_win.exe` | 1233920 | `e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5` | https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_win.exe |
| `vina_split_1.2.7_win.exe` | 687616 | `bc1d4775c1a06749e39818800e9f6d410199c5ab2e0cb136c93f8d53dced8fb3` | https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_split_1.2.7_win.exe |

Para reproducir:

```bash
curl -L -o tools/vina/vina_1.2.7_win.exe https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_win.exe
curl -L -o tools/vina/vina_split_1.2.7_win.exe https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_split_1.2.7_win.exe
sha256sum tools/vina/vina_1.2.7_win.exe tools/vina/vina_split_1.2.7_win.exe
```

Verificar antes de usar en el piloto que los checksums coinciden con los registrados arriba.
