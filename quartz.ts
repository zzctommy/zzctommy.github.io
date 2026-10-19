import { loadQuartzConfig, loadQuartzLayout } from "./quartz/plugins/loader/config-loader"
import { componentRegistry } from "./quartz/components/registry"
import type { ExplorerOptions } from "@quartz-community/explorer"

// npm plugin overrides must use the same package name as quartz.config.yaml.
// The generated .quartz/plugins wrapper uses a different key in this version.
componentRegistry.setOptionOverrides("@quartz-community/explorer", {
  sortFn: (a, b) => {
    if (a.isFolder !== b.isFolder) {
      return a.isFolder ? -1 : 1
    }

    if (a.isFolder) {
      return (a.displayName ?? "").localeCompare(b.displayName ?? "", "zh-CN")
    }

    return (b.slugSegment ?? "").localeCompare(a.slugSegment ?? "", undefined, {
      numeric: true,
    })
  },
} satisfies Partial<ExplorerOptions>)

const config = await loadQuartzConfig()
export default config
export const layout = await loadQuartzLayout()
