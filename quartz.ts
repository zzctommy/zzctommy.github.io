import { loadQuartzConfig, loadQuartzLayout } from "./quartz/plugins/loader/config-loader"
import * as ExternalPlugin from "./.quartz/plugins"

ExternalPlugin.Explorer({
  sortFn: (a, b) => {
    if (a.isFolder !== b.isFolder) {
      return a.isFolder ? -1 : 1
    }

    if (a.isFolder) {
      return (a.displayName ?? "").localeCompare(b.displayName ?? "", "zh-CN")
    }

    return (b.slugSegment ?? "").localeCompare(
      a.slugSegment ?? "",
      undefined,
      { numeric: true },
    )
  },
})

const config = await loadQuartzConfig()
export default config
export const layout = await loadQuartzLayout()